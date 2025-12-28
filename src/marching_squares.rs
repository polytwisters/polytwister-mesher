//! An implementation of Marching Squares algorithm on a cylinder, converting a 2D boolean function
//! R x S^1 -> {0, 1} to a triangular mesh. The function is referred to as "inside," taking
//! arguments u (linear) and theta (circular, from 0 to 2pi) and returning a bool which is true if
//! the point (u, theta) is inside the shape.

use core::f64;
use std::collections::HashMap;
use na::{Point3, Vector3};

use crate::mesh::{self, Face, Mesh, Vertex};

/// Marching Squares mesher for *cylindrical* 2D space with coordinates (u, theta). The grid cells
/// are stored in a 1D vector in a u-major order. For visualization we flatten out the grid and
/// say that u is the vertical direction and theta is the horizontal direction.
pub struct MarchingSquares {
    grid: Grid,
    /// Index u_index * theta_cells + theta_index.
    cells: Vec<Cell>,
}

/// Information about the grid, allowing conversion between discrete and continuous coordinates.
#[derive(Copy, Clone, Debug)]
struct Grid {
    pub u_cells: usize,
    pub theta_cells: usize,
    pub u_min: f64,
    pub u_max: f64,
}

/// One cell in the Marching Squares algorithm.
struct Cell {
    pub shape: u8,
}

/// A unique square in the quad tree. The discrete coordinate of its upper left corner is
/// (u_index, theta_index) / 2^depth and its side length is 2^depth.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
struct MSSquare {
    pub u_index: usize,
    pub theta_index: usize, 
    pub depth: usize,
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

    fn offset(&self, du: isize, dtheta: isize) -> Self {
        MSSquare {
            // Not sure the logic is right here when they become negative
            u_index: (self.u_index as isize + du) as usize,
            theta_index: (self.theta_index as isize + dtheta).rem_euclid(self.theta_cells as isize) as usize,
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
}

impl Grid {
    pub fn u_spacing(&self) -> f64 {
        (self.u_max - self.u_min) / (self.u_cells + 1) as f64
    }

    pub fn theta_spacing(&self) -> f64 {
        f64::consts::TAU / (self.theta_cells + 1) as f64
    }

    /// Convert the index of a U line to its U coordinate. The u_index is the index of the
    /// *boundary* of the cell, not the cell itself.
    pub fn u_index_to_u(&self, u_index: usize) -> f64 {
        let unipolar = u_index as f64 / (self.u_cells + 1) as f64;
        self.u_min + unipolar * (self.u_max - self.u_min)
    }

    /// Convert the index of a theta line to its theta coordinate. The theta_index is the index of
    /// the *boundary* of the cell, not the cell itself.
    pub fn theta_index_to_theta(&self, theta_index: usize) -> f64 {
        (theta_index.rem_euclid(self.theta_cells) as f64 / self.theta_cells as f64) * f64::consts::TAU
    }

    fn vertex_coordinate<F: Fn(f64, f64) -> bool>(&self, vertex: MSPoint, inside: &F) -> (f64, f64) {
        match vertex {
            MSPoint::Corner(square) => (
                self.u_index_to_u(square.u_index),
                self.theta_index_to_theta(square.theta_index)
            ),
            MSPoint::HorizontalEdge(square) => {
                let u_index = square.u_index;
                let theta_index = square.theta_index;
                let t = bisection_search(|t|
                    inside(
                        self.u_index_to_u(u_index),
                        self.theta_index_to_theta(theta_index) + self.theta_spacing() * t
                    ) 
                );
                (
                    self.u_index_to_u(u_index),
                    self.theta_index_to_theta(theta_index) + self.theta_spacing() * t
                )
            },
            MSPoint::VerticalEdge(square) => {
                let u_index = square.u_index;
                let theta_index = square.theta_index;
                let t = bisection_search(|t|
                    inside(
                        self.u_index_to_u(u_index) + self.u_spacing() * t,
                        self.theta_index_to_theta(theta_index)
                    ) 
                );
                (
                    self.u_index_to_u(u_index) + self.u_spacing() * t,
                    self.theta_index_to_theta(theta_index)
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

    fn realize_vertex<F: Fn(f64, f64) -> bool>(&mut self, ms_vertex: MSPoint, grid: &Grid, inside: &F) -> usize {
        if let Some(&mesh_index) = self.ms_vertex_to_mesh_vertex.get(&ms_vertex) {
            return mesh_index;
        }
        let coordinate_2d = grid.vertex_coordinate(ms_vertex, inside);
        let mesh_index = self.mesh_vertices.len();
        self.mesh_vertices.push(Vertex {
            p: Point3::new(
                -coordinate_2d.1.cos(),
                coordinate_2d.1.sin(),
                coordinate_2d.0,
            ),
            n: Vector3::z()
        });
        self.ms_vertex_to_mesh_vertex.insert(ms_vertex, mesh_index);
        mesh_index
    }

    fn into_vertices(self) -> Vec<Vertex> {
        self.mesh_vertices
    }
}

impl MarchingSquares {
    fn new(grid: Grid) -> Self {
        MarchingSquares {
            grid, cells: vec![]
        }
    }

    /// Sample the scalar field and produce a Mesh.
    fn mesh<F: Fn(f64, f64) -> bool>(&self, inside: &F) -> Mesh {
        let mut vertex_index = 0usize;
        let mut realization = MSVertices::new();
        let mut faces = vec![];
        for u in 0..self.grid.u_cells {
            for t in 0..self.grid.theta_cells {
                let triangle_iter = self.mesh_cell(u, t, inside);
                for triangle in triangle_iter {
                    let v1 = realization.realize_vertex(triangle.v1, &self.grid, &inside);
                    let v2 = realization.realize_vertex(triangle.v2, &self.grid, &inside);
                    let v3 = realization.realize_vertex(triangle.v3, &self.grid, &inside);
                    faces.push(Face { v1, v2, v3 });
                }
            }
        }
        Mesh { faces, vertices: realization.into_vertices() }
    }

    /// Sample the shape field at four points and return an appropriate MSTriangle.
    fn mesh_cell<F: Fn(f64, f64) -> bool>(
        &self, 
        u_index: usize,
        theta_index: usize,
        inside: &F
    ) -> impl Iterator<Item = MSTriangle> {
        let iu1 = u_index;
        let iu2 = u_index + 1;
        let it1 = theta_index;
        let it2 = (theta_index + 1).rem_euclid(self.grid.theta_cells);

        let u1 = self.grid.u_index_to_u(u_index);
        let u2 = self.grid.u_index_to_u(u_index + 1);
        let theta1 = self.grid.theta_index_to_theta(theta_index);
        let theta2 = self.grid.theta_index_to_theta(theta_index + 1);

        let square = MSSquare::new_root(u_index, theta_index, self.grid.theta_cells);

        let shape = (
            inside(u1, theta1),
            inside(u1, theta2),
            inside(u2, theta1),
            inside(u2, theta2),
        );

        //    it1     it2
        // iu1 0---1---2
        //     3---x---4
        // iu2 5---6---7
        let v = (
            square.top_left(),
            square.top_middle(),
            square.top_right(),
            square.left_middle(),
            square.right_middle(),
            square.bottom_left(),
            square.bottom_middle(),
            square.bottom_right(),
        );

        match shape {
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
        } .into_iter()
    }
}

#[cfg(test)]
mod test {
    use core::f64;
    use std::path::PathBuf;

    use super::*;

    #[test]
    fn test_ms() {
        let grid = Grid {
            u_cells: 40,
            theta_cells: 40,
            u_min: -2.0,
            u_max: 2.0,
        };
        let ms = MarchingSquares::new(grid);
        let mesh = ms.mesh(&|u, theta|
            u > (theta * 2.0).cos()
        );
        mesh.write_ply_file(&PathBuf::from("out_ms.ply"));
    }
}