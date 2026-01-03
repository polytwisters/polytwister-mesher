#[macro_use]
use nalgebra as na;

use core::{f64, num};
use std::io::{Write};
use std::path::PathBuf;
use std::fs::{File};

use na::{Point, Point3, Transform3, Vector3};
use crate::pipe_section::{PipeSection};

#[derive(Clone, Copy, Debug)]
pub struct Vertex {
    pub p: Point3<f64>,
    pub n: Vector3<f64>,
}

impl Vertex {
    pub fn new(p: Point3<f64>, n: Vector3<f64>) -> Self {
        Vertex { p, n }
    }
}

/// Triangular face with three vertex indices. Vector3<f64> indices start with 0.
#[derive(Clone, Copy, Debug)]
pub struct Face {
    pub v1: usize,
    pub v2: usize,
    pub v3: usize,
}

#[derive(Clone)]
/// A triangular mesh in 3D space. Vertices have normals.
pub struct Mesh {
    pub vertices: Vec<Vertex>,
    pub faces: Vec<Face>,
}


impl Mesh {
    /// Make a new empty mesh.
    pub fn empty() -> Self {
        Mesh { vertices: vec![], faces: vec![] }
    }

    /// A "partial mesh" is one where some vertices may be None. Converting a partial mesh to a mesh
    /// removes the "None" vertices and any faces they are connected to.
    pub fn from_partial(vertices: &Vec<Option<Vertex>>, faces: &Vec<Face>) -> Self {
        let mut new_vertices = vec![];
        let mut new_index = 0usize;
        // Vector of vertex indices whose length is equal to self.vertices.len() such that
        // old_to_new_indices[old_index] is Some(new_vertex_index) if the vertex exists, and None
        // otherwise.
        let mut old_to_new_indices: Vec<Option<usize>> = vec![];
        for option_vertex in vertices.iter() {
            if let Some(vertex) = option_vertex {
                new_vertices.push(vertex.clone());
                old_to_new_indices.push(Some(new_index));
                new_index += 1;
            } else {
                old_to_new_indices.push(None);
            }
        }
        let new_faces = faces.iter().filter_map(|face| {
            let v1_new = old_to_new_indices[face.v1];
            let v2_new = old_to_new_indices[face.v2];
            let v3_new = old_to_new_indices[face.v3];
            if let (Some(v1), Some(v2), Some(v3)) = (v1_new, v2_new, v3_new) {
                Some(Face { v1, v2, v3 })
            } else {
                None
            }
        }).collect::<Vec<_>>();
        Mesh { vertices: new_vertices, faces: new_faces }
    }

    /// Make a spherical mesh.
    pub fn uv_sphere(center: &Point3<f64>, radius: f64, segments: usize, rings: usize) -> Self {
        // Point3<f64> indices: i * segments + j
        // i is the segment index, j is in the ring index.
        let mut vertices: Vec<Vertex> = vec![];
        for i in 0..rings {
            let i_unipolar = (i as f64 + 1.0) / (rings as f64 + 1.0);
            let i_bipolar = i_unipolar * 2.0 - 1.0;
            let elevation = i_bipolar * f64::consts::FRAC_PI_2;
            for j in 0..segments {
                let azimuth = j as f64 / segments as f64 * f64::consts::TAU;
                let normal = Vector3::new(
                    azimuth.cos() * elevation.cos(),
                    azimuth.sin() * elevation.cos(),
                    elevation.sin()
                );
                let point = center + normal * radius;
                let vertex = Vertex::new(point, normal);
                vertices.push(vertex);
            }
        }

        let south_pole_index = vertices.len();
        vertices.push(Vertex::new(
            center + Vector3::new(0.0, 0.0, -radius),
            Vector3::new(0.0, 0.0, -1.0)
        ));
    
        let north_pole_index = vertices.len();
        vertices.push(Vertex::new(
            center + Vector3::new(0.0, 0.0, radius),
            Vector3::new(0.0, 0.0, 1.0)
        ));

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
                v1: v2,
                v2: v1,
                v3: south_pole_index,
            });

            let v3 = (rings - 1) * segments + j;
            let v4 = (rings - 1) * segments + (j + 1) % segments;
            faces.push(Face {
                v1: v3,
                v2: v4,
                v3: north_pole_index,
            });
        }

        Mesh { vertices, faces }
    }

    /// Make a plane parallel to the xy-plane at coordinate z.
    pub fn plane(z: f64, half_length: f64, segments: usize) -> Self {
        Self::partial_plane(|_| true, z, half_length, segments)
    }

    /// Make a plane parallel to the xy-plane at coordinate z, but filter the vertices using a
    /// predicate.
    pub fn partial_plane<F: Fn(&Point3<f64>) -> bool>(predicate: F, z: f64, half_length: f64, segments: usize) -> Self {
        let sign = z.signum();

        // Point3<f64> indices: i * (segments + 1) + j
        let mut vertices: Vec<Option<Vertex>> = vec![];
        for i in 0..=segments {
            let i_unipolar = (i as f64) / (segments as f64);
            let i_bipolar = i_unipolar * 2.0 - 1.0;
            for j in 0..=segments {
                let j_unipolar = (j as f64) / (segments as f64);
                let j_bipolar = j_unipolar * 2.0 - 1.0;
                // Sign flip for negative z to fix handedness of triangles. (It's easier to flip it
                // here than when building the triangles.)
                let x = i_bipolar * half_length * sign;
                let y = j_bipolar * half_length;
                let point = Point3::new(x, y, z);
                let normal = Vector3::new(0.0, 0.0, sign);
                let vertex = Vertex::new(point, normal);
                vertices.push(if predicate(&point) { Some(vertex) } else { None });
            }
        }

        let hop = segments + 1;
        let mut faces: Vec<Face> = vec![];
        for i in 0..segments {
            for j in 0..segments {
                let v1 = i * hop + j;
                let v2 = i * hop + (j + 1) % hop;
                let v3 = (i + 1) * hop + j;
                let v4 = (i + 1) * hop + (j + 1) % hop;

                let v2_exists = matches!(vertices[v2], Some(_));
                let v3_exists = matches!(vertices[v3], Some(_));

                // If all 4 vertices are present, then we triangulate the square
                // with 2 triangles:
                //
                // v1 -- v2
                // | ,--' |
                // v3 -- v4
                //
                // If v2 or v3 are not present, we use the alternate
                // triangulation:
                //
                // v1----v2
                // | `--. |
                // v3 -- v4
                //
                // Mesh::from_partial will delete all triangles with missing
                // vertices.
                //
                // In the above diagrams, the normals face the viewer and the triangles are read off
                // counterclockwise.

                if v2_exists && v3_exists {
                    faces.push(Face {
                        v1: v1,
                        v2: v3,
                        v3: v2,
                    });
                    faces.push(Face {
                        v1: v2,
                        v2: v3,
                        v3: v4,
                    });
                } else {
                    faces.push(Face {
                        v1: v1,
                        v2: v4,
                        v3: v2,
                    });
                    faces.push(Face {
                        v1: v3,
                        v2: v4,
                        v3: v1,
                    });
                }
            }
        }
        Mesh::from_partial(&vertices, &faces)
    }

    pub fn num_vertices(&self) -> usize {
        self.vertices.len()
    }

    pub fn num_faces(&self) -> usize {
        self.faces.len()
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

    /**
     * Given a predicate on vertex locations, return a new Mesh that removes all vertices that do
     * not satisfy that predicate, and any faces that are connected to said vertices.
     */
    pub fn filter_vertices<F: Fn(&Point3<f64>) -> bool>(&self, predicate: F) -> Self {
        let new_vertices = self.vertices.iter().map(|v|
            if (predicate(&v.p)) { Some(v.clone()) } else { None }
        ).collect::<Vec<_>>();
        Mesh::from_partial(&new_vertices, &self.faces)
    }

    pub fn write_obj<W: Write>(&self, buffer: &mut W) -> std::io::Result<()> {
        for vertex in &self.vertices {
            // Blender seems to swap Z and Y for OBJ, so we re-swap them here.
            write!(buffer, "v {} {} {}\n", vertex.p.x, vertex.p.z, vertex.p.y)?;
        }
        for vertex in &self.vertices {
            // Z and Y swapped again here.
            write!(buffer, "vn {} {} {}\n", vertex.n.x, vertex.n.z, vertex.n.y)?;
        }
        for face in &self.faces {
            // OBJ vertex indices start from 1.
            write!(buffer, "f {} {} {}\n", face.v1 + 1, face.v2 + 1, face.v3 + 1)?;
        }
        Ok(())
    }

    pub fn write_ply<W: Write>(&self, buffer: &mut W) -> std::io::Result<()> {
        write!(buffer, "ply\n")?;
        write!(buffer, "format binary_little_endian 1.0\n")?;
        write!(buffer, "element vertex {}\n", self.vertices.len())?;
        write!(buffer, "property float x\n")?;
        write!(buffer, "property float y\n")?;
        write!(buffer, "property float z\n")?;
        write!(buffer, "property float nx\n")?;
        write!(buffer, "property float ny\n")?;
        write!(buffer, "property float nz\n")?;
        write!(buffer, "element face {}\n", self.faces.len())?;
        write!(buffer, "property list uchar int vertex_index\n")?;
        write!(buffer, "end_header\n")?;
        for vertex in &self.vertices {
            buffer.write(&(vertex.p.x as f32).to_le_bytes())?;
            buffer.write(&(vertex.p.y as f32).to_le_bytes())?;
            buffer.write(&(vertex.p.z as f32).to_le_bytes())?;
            buffer.write(&(vertex.n.x as f32).to_le_bytes())?;
            buffer.write(&(vertex.n.y as f32).to_le_bytes())?;
            buffer.write(&(vertex.n.z as f32).to_le_bytes())?;
        }
        for face in &self.faces {
            buffer.write(&[3u8])?;
            // PLY vertices start at 0.
            buffer.write(&(face.v1 as u32).to_le_bytes())?;
            buffer.write(&(face.v2 as u32).to_le_bytes())?;
            buffer.write(&(face.v3 as u32).to_le_bytes())?;
        }
        Ok(())
    }

    pub fn write_ply_file(&self, path: &PathBuf) -> std::io::Result<()> {
        let mut buffer = File::create(path)?;
        self.write_ply(&mut buffer)?;
        Ok(())
    }
}


#[derive(Clone, Copy, Debug)]
pub struct Color {
    pub red: u8,
    pub green: u8,
    pub blue: u8,
}

/// A collection of meshes with an RGB color assigned to each one.
pub struct ColoredMesh {
    pub meshes: Vec<(Mesh, Color)>
}

impl ColoredMesh {
    pub fn write_ply<W: Write>(&self, buffer: &mut W) -> std::io::Result<()> {
        let num_vertices: usize = self.meshes.iter().map(|(mesh, _)| mesh.vertices.len()).sum();
        let num_faces: usize = self.meshes.iter().map(|(mesh, _)| mesh.faces.len()).sum();

        write!(buffer, "ply\n")?;
        write!(buffer, "format binary_little_endian 1.0\n")?;
        write!(buffer, "element vertex {}\n", num_vertices)?;
        write!(buffer, "property float x\n")?;
        write!(buffer, "property float y\n")?;
        write!(buffer, "property float z\n")?;
        write!(buffer, "property float nx\n")?;
        write!(buffer, "property float ny\n")?;
        write!(buffer, "property float nz\n")?;
        write!(buffer, "property uchar red\n")?;
        write!(buffer, "property uchar green\n")?;
        write!(buffer, "property uchar blue\n")?;
        write!(buffer, "element face {}\n", num_faces)?;
        write!(buffer, "property list uchar int vertex_index\n")?;
        write!(buffer, "end_header\n")?;
        for (mesh, color) in &self.meshes {
            for vertex in &mesh.vertices {
                buffer.write(&(vertex.p.x as f32).to_le_bytes())?;
                buffer.write(&(vertex.p.y as f32).to_le_bytes())?;
                buffer.write(&(vertex.p.z as f32).to_le_bytes())?;
                buffer.write(&(vertex.n.x as f32).to_le_bytes())?;
                buffer.write(&(vertex.n.y as f32).to_le_bytes())?;
                buffer.write(&(vertex.n.z as f32).to_le_bytes())?;
                buffer.write(&[color.red])?;
                buffer.write(&[color.green])?;
                buffer.write(&[color.blue])?;
            }
        }
        let mut offset = 0;
        for (mesh, _) in &self.meshes {
            for face in &mesh.faces {
                buffer.write(&[3u8])?;
                // PLY vertices start at 0.
                buffer.write(&((face.v1 + offset) as u32).to_le_bytes())?;
                buffer.write(&((face.v2 + offset) as u32).to_le_bytes())?;
                buffer.write(&((face.v3 + offset) as u32).to_le_bytes())?;
            }
            offset += mesh.vertices.len();
        }
        Ok(())
    }

    pub fn write_ply_file(&self, path: &PathBuf) -> std::io::Result<()> {
        let mut buffer = File::create(path)?;
        self.write_ply(&mut buffer)?;
        Ok(())
    }
}