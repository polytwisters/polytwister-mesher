//! A Marching Squares isosurface mesher specifically for meshing the intersection of a hollow
//! cylinder with an arbitrary set given by an implicit function in R^3. In particular, given the
//! parametric function U : R x S^1 -> R^3 and an indicator function C : R^3 -> {0, 1}, create a
//! triangle mesh that approximates the set of all 3D points x = U(u, theta) such that C(x) = 1.
//! 
//! Although the mesh is in 3D, we do not need the full Marching Cubes algorithm as we are actually
//! meshing an surface that is a subset of an existing 2D surface with an explicit parametrization.

use core::f64;
use std::collections::HashMap;
use na::{Point3, Vector3};

use crate::mesh::{self, Face, Mesh, Vertex};

/// An isosurface comprises two functions: an explicit parametrization that converts surface
/// coordinates (u, theta) into a Vertex with a 3D location and normal, and an indicator function
/// that returns whether the point (u, theta) in surface coordinates is inside the shape.
pub trait Isosurface {
    fn surface_coords_to_point(&self, u: f64, theta: f64) -> Point3<f64>;
    fn surface_coords_to_normal(&self, u: f64, theta: f64) -> Vector3<f64>;
    fn contains_point(&self, p: &Point3<f64>) -> bool;

    fn contains_surface_coord(&self, u: f64, theta: f64) -> bool {
        self.contains_point(&self.surface_coords_to_point(u, theta))
    }
}

/// Marching Squares mesher for *cylindrical* 2D space with coordinates (u, theta). The grid cells
/// are stored in a 1D vector in a u-major order. For visualization we flatten out the grid and
/// say that u is the vertical direction and theta is the horizontal direction.
pub struct MarchingSquares {
    grid: Grid,
    /// Index u_index * theta_cells + theta_index.
    cells: Vec<Node>,
    max_depth: u8,
}

/// Information about the grid, allowing conversion between discrete and continuous coordinates.
#[derive(Copy, Clone, Debug)]
pub struct Grid {
    pub u_cells: usize,
    pub theta_cells: usize,
    pub u_min: f64,
    pub u_max: f64,
}

/// A node in the quadtree. It knows its gometric location as an MSSquare.
struct Node {
    square: MSSquare,
    subtree: Subtree,
}

/// The subtree associated with a node.
enum Subtree {
    /// The node has four children.
    Parent(Box<(Node, Node, Node, Node)>),
    /// A square at the lowest level of the quadtree which contains the boundary of the shape. The
    /// four bools indicate which of the top left, top right, bottom left, and bottom right corners
    /// are in the shape.
    Leaf((bool, bool, bool, bool)),
    /// A square entirely inside the shape.
    Full,
    /// A square entirely outside the shape.
    Empty
}

/// A unique square in the quad tree. The discrete coordinate of its upper left corner is
/// (u_index, theta_index) / 2^depth and its side length is 2^depth.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
struct MSSquare {
    pub u_index: usize,
    pub theta_index: usize, 
    pub depth: u8,
    // The MSSquare does need to know the number of radial segments of the cylinder.
    pub theta_cells: usize,
}

/// A discrete vertex location used during marching squares. MSPoint::Corner is a vertex located
/// at the corner of a square. HorizontalEdge and VerticalEdge are located on edges of squares.
/// In the latter two cases, the floating-point coordinates are not known yet, and are later derived
/// using interpolation.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
enum MSPoint {
    Corner(MSSquare),
    HorizontalEdge(MSSquare),
    VerticalEdge(MSSquare),
}

struct MSTriangle {
    pub v1: MSPoint,
    pub v2: MSPoint,
    pub v3: MSPoint,
}

/// A collection of marching squares vertices, maintaining their "MS coordinates" (discrete coords
/// with information on whether vertices are on corners/edges) and Euclidean coordinates and
/// indices in the final Mesh object to be produced.
struct MSVertices {
    ms_vertex_to_mesh_vertex: HashMap<MSPoint, usize>,
    mesh_vertices: Vec<Vertex>,
}

/// Given a monotonic function f: [0, 1] -> bool, use bisection search to find an x in [0, 1] so
/// f(x) is right on the cusp of the switch from "true" to "false" or vice versa.
fn bisection_search<F: Fn(f64) -> bool>(f: F) -> f64 {
    let mut x_min = 0.0; 
    let mut x_max = 1.0;
    if f(x_min) == f(x_max) {
        return 1.0;
    }
    // True if the function is ramping from false to true.
    let upward = f(x_max);
    let mut x;
    for i in 0..5 {
        x = (x_min + x_max) / 2.0;
        // This flips the if/else statements if upward is false.
        if f(x) == upward {
            x_max = x;
        } else {
            x_min = x;
        }
    }
    (x_min + x_max) / 2.0
}

impl MSSquare {
    fn new_root(u_index: usize, theta_index: usize, theta_cells: usize) -> Self {
        MSSquare { u_index, theta_index, theta_cells, depth: 0 }
    }

    /// Offset this MSSquare to another MSSquare of the same depth. Only positive offsets are
    /// allowed.
    fn offset(&self, du: usize, dtheta: usize) -> Self {
        let theta_index = ((self.theta_index + dtheta) as u64).rem_euclid(
            self.theta_cells as u64 * Grid::inv_depth_scale(self.depth)
        );
        MSSquare {
            u_index: (self.u_index + du) as usize,
            theta_index: theta_index as usize,
            depth: self.depth,
            theta_cells: self.theta_cells,
        }
    }

    fn top_left(&self) -> MSPoint {
        MSPoint::Corner(self.clone())
    }

    fn top_right(&self) -> MSPoint {
        MSPoint::Corner(self.offset(0, 1))
    }

    fn bottom_left(&self) -> MSPoint {
        MSPoint::Corner(self.offset(1, 0))
    }

    fn bottom_right(&self) -> MSPoint {
        MSPoint::Corner(self.offset(1, 1))
    }

    fn top_middle(&self) -> MSPoint {
        MSPoint::HorizontalEdge(self.clone())
    }

    fn left_middle(&self) -> MSPoint {
        MSPoint::VerticalEdge(self.clone())
    }

    fn bottom_middle(&self) -> MSPoint {
        MSPoint::HorizontalEdge(self.offset(1, 0))
    }

    fn right_middle(&self) -> MSPoint {
        MSPoint::VerticalEdge(self.offset(0, 1))
    }

    /// Given the "inside-ness" of the four corners of this square (true if inside the shape and
    /// false otherwise), triangulate the shape in this square.
    fn triangulate(
        &self,
        corners: (bool, bool, bool, bool)
    ) -> Vec<MSTriangle> {
        //    it1     it2
        // iu1 0---1---2
        //     3---x---4
        // iu2 5---6---7
        let v = (
            self.top_left(),
            self.top_middle(),
            self.top_right(),
            self.left_middle(),
            self.right_middle(),
            self.bottom_left(),
            self.bottom_middle(),
            self.bottom_right(),
        );

        match corners {
            // Empty
            (false, false, false, false) => vec![],
    
            // Full
            (true, true, true, true) => vec![
                MSTriangle { v1: v.0, v2: v.5, v3: v.2 },
                MSTriangle { v1: v.2, v2: v.5, v3: v.7 },
            ],

            // Single corner
            (true, false, false, false) => vec![
                MSTriangle { v1: v.0, v2: v.3, v3: v.1 },
            ],
            (false, true, false, false) => vec![
                MSTriangle { v1: v.2, v2: v.1, v3: v.4 },
            ],
            (false, false, true, false) => vec![
                MSTriangle { v1: v.3, v2: v.5, v3: v.6 },
            ],
            (false, false, false, true) => vec![
                MSTriangle { v1: v.4, v2: v.6, v3: v.7 },
            ],

            // Opposite, diagonal corners
            (true, false, false, true) => vec![
                MSTriangle { v1: v.0, v2: v.3, v3: v.1 },
                MSTriangle { v1: v.4, v2: v.6, v3: v.7 },
            ],
            (false, true, true, false) => vec![
                MSTriangle { v1: v.2, v2: v.1, v3: v.4 },
                MSTriangle { v1: v.3, v2: v.5, v3: v.6 },
            ],

            // Two adjacent corners
            (true, true, false, false) => vec![
                MSTriangle { v1: v.0, v2: v.3, v3: v.2 },
                MSTriangle { v1: v.3, v2: v.4, v3: v.2 },
            ],
            (true, false, true, false) => vec![
                MSTriangle { v1: v.0, v2: v.5, v3: v.1 },
                MSTriangle { v1: v.1, v2: v.5, v3: v.6 },
            ],
            (false, false, true, true) => vec![
                MSTriangle { v1: v.3, v2: v.5, v3: v.4 },
                MSTriangle { v1: v.5, v2: v.7, v3: v.4 },
            ],
            (false, true, false, true) => vec![
                MSTriangle { v1: v.1, v2: v.6, v3: v.2 },
                MSTriangle { v1: v.6, v2: v.7, v3: v.2 },
            ],

            // Three corners
            (true, true, true, false) => vec![
                MSTriangle { v1: v.0, v2: v.4, v3: v.2 },
                MSTriangle { v1: v.0, v2: v.6, v3: v.4 },
                MSTriangle { v1: v.0, v2: v.5, v3: v.6 },
            ],
            (true, false, true, true) => vec![
                MSTriangle { v1: v.0, v2: v.5, v3: v.1 },
                MSTriangle { v1: v.1, v2: v.5, v3: v.4 },
                MSTriangle { v1: v.4, v2: v.5, v3: v.7 },
            ],
            (true, true, false, true) => vec![
                MSTriangle { v1: v.2, v2: v.0, v3: v.3 },
                MSTriangle { v1: v.2, v2: v.3, v3: v.6 },
                MSTriangle { v1: v.2, v2: v.6, v3: v.7 },
            ],
            (false, true, true, true) => vec![
                MSTriangle { v1: v.7, v2: v.2, v3: v.1 },
                MSTriangle { v1: v.7, v2: v.1, v3: v.3 },
                MSTriangle { v1: v.7, v2: v.3, v3: v.5 },
            ],
            _ => panic!()
        }
    }

    fn subdivide(&self) -> (MSSquare, MSSquare, MSSquare, MSSquare) {
        let top_left = MSSquare {
            u_index: self.u_index * 2,
            theta_index: self.theta_index * 2,
            depth: self.depth + 1,
            theta_cells: self.theta_cells,
        };
        (
            top_left,
            top_left.offset(0, 1),
            top_left.offset(1, 0),
            top_left.offset(1, 1),
        )
    }
}

impl Node {
    fn triangulate(&self) -> Vec<MSTriangle> {
        match &self.subtree {
            Subtree::Empty => vec![],
            Subtree::Full => self.square.triangulate((true, true, true, true)),
            Subtree::Leaf(corners) => self.square.triangulate(corners.clone()),
            Subtree::Parent(children) => {
                let mut result = vec![];
                result.extend(children.0.triangulate());
                result.extend(children.1.triangulate());
                result.extend(children.2.triangulate());
                result.extend(children.3.triangulate());
                result
            }
        }
    }
}

impl Grid {
    fn inv_depth_scale(depth: u8) -> u64 {
        1 << (depth as u64)
    }

    fn depth_scale(depth: u8) -> f64 {
        1.0 / Grid::inv_depth_scale(depth) as f64
    }

    pub fn u_spacing(&self) -> f64 {
        (self.u_max - self.u_min) / (self.u_cells + 1) as f64
    }

    pub fn theta_spacing(&self) -> f64 {
        f64::consts::TAU / (self.theta_cells + 1) as f64
    }

    /// Convert the index of a U line to its U coordinate. The u_index is the index of the
    /// *boundary* of the cell, not the cell itself.
    pub fn u_index_to_u(&self, u_index: usize, depth: u8) -> f64 {
        self.u_index_to_u_continuous(u_index as f64, depth)
    }

    /// Convert the index of a theta line to its theta coordinate. The theta_index is the index of
    /// the *boundary* of the cell, not the cell itself.
    pub fn theta_index_to_theta(&self, theta_index: usize, depth: u8) -> f64 {
        self.theta_index_to_theta_continuous(theta_index as f64, depth)
    }

    pub fn u_index_to_u_continuous(&self, u_index: f64, depth: u8) -> f64 {
        let unipolar = u_index / (self.u_cells + 1) as f64 * Grid::depth_scale(depth);
        self.u_min + unipolar * (self.u_max - self.u_min)
    }

    pub fn theta_index_to_theta_continuous(&self, theta_index: f64, depth: u8) -> f64 {
        let tmp = theta_index / self.theta_cells as f64;
        (tmp * Grid::depth_scale(depth)).rem_euclid(1.0) * f64::consts::TAU
    }

    fn vertex_coordinate(&self, vertex: MSPoint, isosurface: &impl Isosurface) -> (f64, f64) {
        match vertex {
            MSPoint::Corner(square) => (
                self.u_index_to_u(square.u_index, square.depth),
                self.theta_index_to_theta(square.theta_index, square.depth)
            ),
            MSPoint::HorizontalEdge(square) => {
                let u_index = square.u_index;
                let theta_index = square.theta_index;
                let depth = square.depth;
                // TODO fix code dupe ewww
                let t = bisection_search(|t|
                    isosurface.contains_surface_coord(
                        self.u_index_to_u(u_index, depth),
                        self.theta_index_to_theta_continuous(theta_index as f64 + t, depth)
                    ) 
                );
                (
                    self.u_index_to_u(u_index, depth),
                    self.theta_index_to_theta_continuous(theta_index as f64 + t, depth)
                )
            },
            MSPoint::VerticalEdge(square) => {
                let u_index = square.u_index;
                let theta_index = square.theta_index;
                let depth = square.depth;
                let t = bisection_search(|t|
                    isosurface.contains_surface_coord(
                        self.u_index_to_u_continuous(u_index as f64 + t, depth),
                        self.theta_index_to_theta(theta_index, depth)
                    ) 
                );
                (
                    self.u_index_to_u_continuous(u_index as f64 + t, depth),
                    self.theta_index_to_theta(theta_index, depth)
                )
            },
        }
    }
}

impl MSVertices {
    fn new() -> Self {
        MSVertices {
            ms_vertex_to_mesh_vertex: HashMap::new(),
            mesh_vertices: vec![],
        }
    }

    fn realize_vertex(
        &mut self,
        ms_vertex: MSPoint,
        grid: &Grid,
        isosurface: &impl Isosurface
    ) -> usize {
        if let Some(&mesh_index) = self.ms_vertex_to_mesh_vertex.get(&ms_vertex) {
            return mesh_index;
        }
        let coordinate_2d = grid.vertex_coordinate(ms_vertex, isosurface);
        let mesh_index = self.mesh_vertices.len();
        self.mesh_vertices.push(Vertex {
            p: isosurface.surface_coords_to_point(coordinate_2d.0, coordinate_2d.1),
            n: isosurface.surface_coords_to_normal(coordinate_2d.0, coordinate_2d.1),
        });
        self.ms_vertex_to_mesh_vertex.insert(ms_vertex, mesh_index);
        mesh_index
    }

    fn into_vertices(self) -> Vec<Vertex> {
        self.mesh_vertices
    }
}

impl MarchingSquares {
    fn new(grid: Grid, max_depth: u8) -> Self {
        MarchingSquares {
            grid,
            cells: Vec::with_capacity(grid.theta_cells * grid.u_cells),
            max_depth,
        }
    }

    fn make_node(
        &self, 
        square: MSSquare,
        isosurface: &impl Isosurface
    ) -> Node {
        // TODO: Move this logic to MSSquare.
        let iu1 = square.u_index;
        let iu2 = square.u_index + 1;
        let it1 = square.theta_index;
        let it2 = (square.theta_index + 1).rem_euclid(
            self.grid.theta_cells * Grid::inv_depth_scale(square.depth) as usize
        );

        let u1 = self.grid.u_index_to_u(iu1, square.depth);
        let u2 = self.grid.u_index_to_u(iu2, square.depth);
        let theta1 = self.grid.theta_index_to_theta(it1, square.depth);
        let theta2 = self.grid.theta_index_to_theta(it2, square.depth);

        let corners = (
            isosurface.contains_surface_coord(u1, theta1),
            isosurface.contains_surface_coord(u1, theta2),
            isosurface.contains_surface_coord(u2, theta1),
            isosurface.contains_surface_coord(u2, theta2),
        );

        let subtree = match corners {
            (true, true, true, true) => Subtree::Full,
            (false, false, false, false) => Subtree::Empty,
            _ => {
                if square.depth >= self.max_depth {
                    Subtree::Leaf(corners)
                } else {
                    let subsquares = square.subdivide();
                    Subtree::Parent(
                        Box::new((
                            self.make_node(subsquares.0, isosurface),
                            self.make_node(subsquares.1, isosurface),
                            self.make_node(subsquares.2, isosurface),
                            self.make_node(subsquares.3, isosurface),
                        ))
                    )
                }
            }
        };

        Node {
            square,
            subtree
        }
    }

    /// Sample the isosurface and produce a Mesh.
    fn mesh(&mut self, isosurface: &impl Isosurface) -> Mesh {
        for u_index in 0..self.grid.u_cells {
            for t_index in 0..self.grid.theta_cells {
                let square = MSSquare::new_root(u_index, t_index, self.grid.theta_cells);
                self.cells.push(self.make_node(square, isosurface));
            }
        }

        let mut realization = MSVertices::new();
        let mut faces = vec![];
        for cell in self.cells.iter() {
            for triangle in cell.triangulate() {
                let v1 = realization.realize_vertex(triangle.v1, &self.grid, isosurface);
                let v2 = realization.realize_vertex(triangle.v2, &self.grid, isosurface);
                let v3 = realization.realize_vertex(triangle.v3, &self.grid, isosurface);
                faces.push(Face { v1, v2, v3 });
            }
        }
        Mesh { faces, vertices: realization.into_vertices() }
    }
}

pub fn meshify(isosurface: &impl Isosurface, grid: &Grid, max_depth: u8) -> Mesh {
    let mut ms = MarchingSquares::new(*grid, max_depth);
    ms.mesh(isosurface)
}

#[cfg(test)]
mod test {
    use core::f64;
    use std::path::PathBuf;

    use na::Point2;

    use super::*;

    struct ExampleIsosurface;

    impl Isosurface for ExampleIsosurface {
        fn surface_coords_to_point(&self, u: f64, theta: f64) -> Point3<f64> {
            Point3::new(theta.cos(), theta.sin(), u)
        }
        fn surface_coords_to_normal(&self, u: f64, theta: f64) -> Vector3<f64> {
            Vector3::new(theta.cos(), theta.sin(), 0.0)
        }
        fn contains_point(&self, p: &Point3<f64>) -> bool {
            p.z < p.x.sin()
        }
    }

    #[test]
    fn test_ms() {
        let grid = Grid {
            u_cells: 40,
            theta_cells: 40,
            u_min: -2.0,
            u_max: 2.0,
        };
        let mesh = meshify(&ExampleIsosurface { }, &grid, 0);
        mesh.write_ply_file(&PathBuf::from("out_ms.ply"));
    }
}