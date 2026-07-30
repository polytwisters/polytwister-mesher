//! A Marching Squares mesher for 2D manifolds embedded in 3D space.
//! 
//! In its basic usage, this mesher takes a 2D surface given as an explicit parametrization
//! U(u, v) : R^2 -> R^3 and an implicit indicator function C : R^3 -> {0, 1}, creating a triangle
//! mesh in 3D space that approximates the 2D surface formed by all 3D points x = U(u, v) such that
//! C(x) = 1. In addition, either u or v can be a circular variable, so the surface U can be a
//! topological cylinder or torus.
//! 
//! This is not the full 3D Marching Cubes algorithm, but rather the simpler 2D version of it,
//! which happens to be on a surface that is in 3D space and possibly curved.

use core::f32;
use std::collections::HashMap;
use na::{Point3, Vector3};

use crate::mesh::{self, Face, Mesh, Vertex};
use crate::utils::bisection_search;

mod grid;
pub use grid::{Grid, GridAxis};

/// An isosurface comprises two functions: an explicit parametrization that converts surface
/// coordinates (u, v) into a Vertex with a 3D location and normal, and an indicator function
/// that returns whether the point (u, v) in surface coordinates is inside the shape.
pub trait Isosurface {
    fn surface_coords_to_point(&self, u: f32, v: f32) -> Point3<f32>;
    fn surface_coords_to_normal(&self, u: f32, v: f32) -> Vector3<f32>;
    fn contains_point(&self, p: &Point3<f32>) -> bool;

    fn contains_surface_coord(&self, u: f32, v: f32) -> bool {
        self.contains_point(&self.surface_coords_to_point(u, v))
    }
}

/// Marching Squares mesher for 2D space with coordinates (u, v). The grid cells
/// are stored in a 1D vector in a u-major order. For visualization we flatten out the grid and
/// say that u is the vertical direction and v is the horizontal direction.
pub struct MarchingSquares {
    grid: Grid,
    cells: Vec<Cell>,
}

/// A cell in Marching Squares whose corners have been determined to be inside or outside the
/// implicit surface.
struct Cell {
    square: Square,
    corners: (bool, bool, bool, bool),
}

/// A unit-size square in Marching Squares.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
struct Square {
    /// u index of the top left corner.
    pub ui: usize,
    /// v index of the top left corner.
    pub vi: usize,
}

/// A point on the grid which is either a corner or on the edge of a Square. If it's on an edge,
/// we do not know at this time where it is on the edge -- that will be computed later.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
enum MSPoint {
    Corner(Square),
    HorizontalEdge(Square),
    VerticalEdge(Square),
}

/// A triangle produced by Marching Squares.
struct MSTriangle {
    pub v1: MSPoint,
    pub v2: MSPoint,
    pub v3: MSPoint,
}

/// A collection of marching triangle vertices, with a mapping from MSPoints to vertex indices in the
/// final Mesh object.
struct MSVertices {
    ms_vertex_to_mesh_vertex: HashMap<MSPoint, usize>,
    mesh_vertices: Vec<Vertex>,
}

impl Square {
    fn new(ui: usize, vi: usize) -> Self {
        Square { ui, vi }
    }

    /// Given a Grid, produce a new square which wraps the indices as necessary for circular axes.
    fn wrap(&self, grid: &Grid) -> Self {
        Square {
            ui: grid.u_axis.wrap(self.ui),
            vi: grid.v_axis.wrap(self.vi),
        }
    }

    /// Offset this Square to another Square. Only positive offsets are allowed. The units of the
    /// offset are such that offsets of (1, 0) and (0, 1) move to vertically and horizontally
    /// adjacent squares.
    fn offset(&self, du: usize, dv: usize) -> Self {
        Square {
            ui: self.ui + du,
            vi: self.vi + dv,
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
}

impl MSPoint {
    /// Given a Grid, produce a new MSPoint which wraps the indices as necessary for circular axes.
    fn wrap(&self, grid: &Grid) -> Self {
        match self {
            Self::Corner(square) => Self::Corner(square.wrap(grid)),
            Self::HorizontalEdge(square) => Self::HorizontalEdge(square.wrap(grid)),
            Self::VerticalEdge(square) => Self::VerticalEdge(square.wrap(grid)),
        }
    }
}

impl Cell {
    fn triangulate(&self) -> Vec<MSTriangle> {
        self.square.triangulate(self.corners.clone())
    }
}

impl Grid {
    fn vertex_coordinate(&self, vertex: MSPoint, isosurface: &impl Isosurface) -> (f32, f32) {
        match vertex {
            MSPoint::Corner(square) => (
                self.ui_to_u(square.ui as f32),
                self.vi_to_v(square.vi as f32)
            ),
            MSPoint::HorizontalEdge(square) => {
                let ui = square.ui as f32;
                let vi = square.vi as f32;
                // TODO fix code dupe ewww
                let t = bisection_search(|t|
                    isosurface.contains_surface_coord(
                        self.ui_to_u(ui),
                        self.vi_to_v(vi + t)
                    ) 
                );
                (
                    self.ui_to_u(ui),
                    self.vi_to_v(vi + t)
                )
            },
            MSPoint::VerticalEdge(square) => {
                let ui = square.ui as f32;
                let vi = square.vi as f32;
                // TODO fix code dupe ewww
                let t = bisection_search(|t|
                    isosurface.contains_surface_coord(
                        self.ui_to_u(ui + t),
                        self.vi_to_v(vi)
                    ) 
                );
                (
                    self.ui_to_u(ui + t),
                    self.vi_to_v(vi)
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
        // Wrap the incoming vertex in case of cylindrical topology.
        let ms_vertex_normalized = ms_vertex.wrap(grid);
        // Try to retrieve a cached vertex if possible.
        if let Some(&mesh_index) = self.ms_vertex_to_mesh_vertex.get(&ms_vertex_normalized) {
            return mesh_index;
        }
        // Create a new Vertex, using the isosurface methods to compute its location and normal.
        let coordinate_2d = grid.vertex_coordinate(ms_vertex_normalized, isosurface);
        let mesh_index = self.mesh_vertices.len();
        self.mesh_vertices.push(Vertex {
            p: isosurface.surface_coords_to_point(coordinate_2d.0, coordinate_2d.1),
            n: isosurface.surface_coords_to_normal(coordinate_2d.0, coordinate_2d.1),
        });
        // Cache its index and return.
        self.ms_vertex_to_mesh_vertex.insert(ms_vertex_normalized, mesh_index);
        mesh_index
    }

    fn into_vertices(self) -> Vec<Vertex> {
        self.mesh_vertices
    }
}

impl MarchingSquares {
    fn new(grid: Grid) -> Self {
        let capacity = grid.u_axis.num_segments() * grid.v_axis.num_segments();
        MarchingSquares {
            grid,
            cells: Vec::with_capacity(capacity),
        }
    }

    fn make_node(
        &self, 
        square: Square,
        isosurface: &impl Isosurface
    ) -> Cell {
        // TODO: Move this logic to MSSquare.
        let ui1 = square.ui;
        let ui2 = square.ui + 1;
        let vi1 = square.vi;
        let vi2 = square.vi + 1;

        let u1 = self.grid.ui_to_u(ui1 as f32);
        let u2 = self.grid.ui_to_u(ui2 as f32);
        let v1 = self.grid.vi_to_v(vi1 as f32);
        let v2 = self.grid.vi_to_v(vi2 as f32);

        let corners = (
            isosurface.contains_surface_coord(u1, v1),
            isosurface.contains_surface_coord(u1, v2),
            isosurface.contains_surface_coord(u2, v1),
            isosurface.contains_surface_coord(u2, v2),
        );

        Cell {
            square,
            corners
        }
    }

    /// Sample the isosurface and produce a Mesh.
    fn mesh(&mut self, isosurface: &impl Isosurface) -> Mesh {
        for ui in 0..self.grid.u_axis.num_segments() {
            for vi in 0..self.grid.v_axis.num_segments() {
                let square = Square::new(ui, vi);
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

pub fn meshify(isosurface: &impl Isosurface, grid: Grid) -> Mesh {
    let mut ms = MarchingSquares::new(grid);
    ms.mesh(isosurface)
}

#[cfg(test)]
mod test {
    use core::f32;
    use std::path::PathBuf;

    use na::Point2;

    use crate::mesh::MeshLike;
    use super::*;

    struct ExampleIsosurface;

    impl Isosurface for ExampleIsosurface {
        fn surface_coords_to_point(&self, u: f32, v: f32) -> Point3<f32> {
            Point3::new(-u.cos(), u.sin(), v)
        }
        fn surface_coords_to_normal(&self, u: f32, v: f32) -> Vector3<f32> {
            Vector3::new(-u.cos(), u.sin(), 0.0)
        }
        fn contains_point(&self, p: &Point3<f32>) -> bool {
            p.z < (p.x * 8.0).sin() * 0.5
        }
    }

    #[test]
    fn test_cylinder() {
        let grid = Grid {
            u_axis: GridAxis::uniform_circular(20, f32::consts::TAU),
            v_axis: GridAxis::uniform_linear(30, -2.0, 2.0),
        };
        let mesh = meshify(&ExampleIsosurface { }, grid);
        mesh.write_ply_file(&PathBuf::from("cylinder.ply"));
    }

    struct ExampleIsosurface2;

    impl Isosurface for ExampleIsosurface2 {
        fn surface_coords_to_point(&self, u: f32, v: f32) -> Point3<f32> {
            Point3::new(u, v, 0.0)
        }
        fn surface_coords_to_normal(&self, u: f32, v: f32) -> Vector3<f32> {
            Vector3::new(0.0, 0.0, 1.0)
        }
        fn contains_point(&self, p: &Point3<f32>) -> bool {
            p.x.hypot(p.y) < 1.0
        }
    }

    #[test]
    fn test_plane() {
        let grid = Grid {
            u_axis: GridAxis::uniform_linear(20, -2.0, 2.0),
            v_axis: GridAxis::uniform_linear(30, -2.0, 2.0),
        };
        let mesh = meshify(&ExampleIsosurface2 { }, grid);
        mesh.write_ply_file(&PathBuf::from("plane.ply"));
    }

    struct ExampleIsosurface3;

    impl Isosurface for ExampleIsosurface3 {
        fn surface_coords_to_point(&self, u: f32, v: f32) -> Point3<f32> {
            Point3::new(u.cos(), u.sin(), v)
        }
        fn surface_coords_to_normal(&self, u: f32, v: f32) -> Vector3<f32> {
            Vector3::new(u.cos(), u.sin(), 0.0)
        }
        fn contains_point(&self, p: &Point3<f32>) -> bool {
            true
        }
    }

    #[test]
    fn test_triangular_prism() {
        let grid = Grid {
            u_axis: GridAxis::uniform_circular(3, f32::consts::TAU),
            v_axis: GridAxis::uniform_linear(2, -2.0, 2.0),
        };
        let mesh = meshify(&ExampleIsosurface3 { }, grid);
        mesh.write_ply_file(&PathBuf::from("prism.ply"));
        // This should generate a triangular prism (with open caps) which has 6 vertices. If it has
        // 8 vertices, then something is wrong with how the circular axis is treated.
        assert_eq!(mesh.num_vertices(), 6);
    }
}