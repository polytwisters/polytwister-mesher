#[macro_use]
use nalgebra as na;

use core::f64;
use std::io::{Write};

use na::{Point, Point3, Transform3, Vector3};
use crate::pipe_section::{PipeSection};

/**
 * Triangular face with three vertex indices. Vector3<f64> indices start with 0.
 */
#[derive(Clone, Copy, Debug)]
pub struct Face {
    pub v1: usize,
    pub v2: usize,
    pub v3: usize,
}

pub struct Mesh {
    pub vertices: Vec<Point3<f64>>,
    pub faces: Vec<Face>,
}


impl Mesh {
    pub fn empty() -> Self {
        Mesh { vertices: vec![], faces: vec![] }
    }

    pub fn uv_sphere(center: &Point3<f64>, radius: f64, segments: usize, rings: usize) -> Self {
        // Point3<f64> indices: i * segments + j
        // i is the segment index, j is in the ring index.
        let mut vertices: Vec<Point3<f64>> = vec![];
        for i in 0..rings {
            let i_unipolar = (i as f64 + 1.0) / (rings as f64 + 1.0);
            let i_bipolar = i_unipolar * 2.0 - 1.0;
            let elevation = i_bipolar * f64::consts::FRAC_PI_2;
            for j in 0..segments {
                let azimuth = j as f64 / segments as f64 * f64::consts::TAU;
                vertices.push(center + Vector3::new(
                    azimuth.cos() * elevation.cos(),
                    azimuth.sin() * elevation.cos(),
                    elevation.sin()
                ) * radius);
            }
        }

        let south_pole_index = vertices.len();
        vertices.push(center + Vector3::new(0.0, 0.0, -1.0));
    
        let north_pole_index = vertices.len();
        vertices.push(center + Vector3::new(0.0, 0.0, 1.0));

        // Connect everything except poles with cylinder topology.
        let mut faces: Vec<Face> = vec![];
        for i in 0..rings - 1 {
            for j in 0..segments {
                let v1 = i * segments + j;
                let v2 = i * segments + (j + 1) % segments;
                let v3 = (i + 1) * segments + j;
                let v4 = (i + 1) * segments + (j + 1) % segments;

                // ^ north pole
                //
                // v3 -- v4
                // | '--. |
                // v1 -- v2
                //
                // v south pole
                faces.push(Face {
                    v1: v1,
                    v2: v2,
                    v3: v3,
                });
                faces.push(Face {
                    v1: v2,
                    v2: v4,
                    v3: v3,
                });
            }
        }

        for j in 0..segments {
            // north pole
            //  /    \
            // v3 -- v4
            //
            // v1 -- v2
            //  \    /
            // south pole

            let v1 = j;
            let v2 = (j + 1) % segments;
            faces.push(Face {
                v1: v1,
                v2: v2,
                v3: south_pole_index,
            });

            let v3 = (rings - 1) * segments + j;
            let v4 = (rings - 1) * segments + (j + 1) % segments;
            faces.push(Face {
                v1: v4,
                v2: v3,
                v3: north_pole_index,
            });
        }

        Mesh { vertices, faces }
    }

    /**
     * Consume this mesh and transform it into a new one using the given 3D transform.
     */
    fn transform(self, transform: &Transform3<f64>) -> Mesh {
        let new_vertices = self.vertices.into_iter().map(|vertex|
            transform.transform_point(&vertex)
        ).collect::<Vec<_>>();
        Mesh { vertices: new_vertices, faces: self.faces }
    }

    /**
     * Make a plane parallel to the xy-plane at coordinate z.
     */
    pub fn plane(z: f64, half_length: f64, segments: usize) -> Self {
        // Point3<f64> indices: i * segments + j
        let mut vertices: Vec<Point3<f64>> = vec![];
        for i in 0..=segments {
            let i_unipolar = (i as f64) / (segments as f64);
            let i_bipolar = i_unipolar * 2.0 - 1.0;
            for j in 0..=segments {
                let j_unipolar = (j as f64) / (segments as f64);
                let j_bipolar = j_unipolar * 2.0 - 1.0;
                let x = i_bipolar * half_length;
                let y = j_bipolar * half_length;
                vertices.push(Point3::new(x, y, z));
            }
        }

        let mut faces: Vec<Face> = vec![];
        for i in 0..segments {
            for j in 0..segments {
                let v1 = i * segments + j;
                let v2 = i * segments + (j + 1);
                let v3 = (i + 1) * segments + j;
                let v4 = (i + 1) * segments + (j + 1);

                // v1 -- v2
                // | ,--' |
                // v3 -- v4
                faces.push(Face {
                    v1: v1,
                    v2: v2,
                    v3: v3,
                });
                faces.push(Face {
                    v1: v2,
                    v2: v4,
                    v3: v3,
                });
            }
        }
        Mesh { vertices, faces }
    }

    fn offset_indices(faces: Vec<Face>, offset: usize) -> Vec<Face> {
        faces.into_iter().map(|face| {
            Face {
                v1: face.v1 + offset,
                v2: face.v2 + offset,
                v3: face.v3 + offset,
            }
        }).collect::<Vec<_>>()
    }

    pub fn merge(meshes: Vec<Self>) -> Self {
        let mut vertices = vec![];
        let mut faces = vec![];
        let mut offset = 0usize;
        for mesh in meshes {
            let num_vertices = mesh.vertices.len();
            faces.extend(Mesh::offset_indices(mesh.faces, offset));
            vertices.extend(mesh.vertices);
            offset += num_vertices;
        }
        Mesh { vertices, faces }
    }

    pub fn write_obj<W: Write>(&self, buffer: &mut W) -> std::io::Result<()> {
        for vertex in &self.vertices {
            write!(buffer, "v {} {} {}\n", vertex.x, vertex.y, vertex.z)?;
        }
        for face in &self.faces {
            // OBJ vertex indices start from 1.
            write!(buffer, "f {} {} {}\n", face.v1 + 1, face.v2 + 1, face.v3 + 1)?;
        }
        Ok(())
    }

    /**
     * Given a predicate on vertex locations, return a new Mesh that removes all vertices that do
     * not satisfy that predicate, and any faces that are connected to said vertices.
     */
    pub fn filter_vertices<F: Fn(&Point3<f64>) -> bool>(&self, predicate: F) -> Self {
        let mut vertices = vec![];
        let mut new_index = 0usize;
        // Vector of vertex indices whose length is equal to self.vertices.len() such that
        // old_to_new_indices[old_index] is Some(new_vertex_index) if the vertex is kept, and None
        // otherwise.
        let mut old_to_new_indices: Vec<Option<usize>> = vec![];
        for vertex in self.vertices.iter() {
            if predicate(&vertex) {
                vertices.push(vertex.clone());
                old_to_new_indices.push(Some(new_index));
                new_index += 1;
            } else {
                old_to_new_indices.push(None);
            }
        }
        let faces = self.faces.iter().filter_map(|face| {
            let v1_new = old_to_new_indices[face.v1];
            let v2_new = old_to_new_indices[face.v2];
            let v3_new = old_to_new_indices[face.v3];
            if let (Some(v1), Some(v2), Some(v3)) = (v1_new, v2_new, v3_new) {
                Some(Face { v1, v2, v3 })
            } else {
                None
            }
        }).collect::<Vec<_>>();
        Mesh { vertices, faces }
    }
}