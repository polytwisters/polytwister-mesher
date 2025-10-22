use std::{fs::File, io::Write};

pub struct Pipe {
    pub a: f32,
    pub b: f32,
    pub c: f32,
}

/**
 * Mesh vertex, location given in Cartesian coordinates.
 */
pub struct Vertex {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

/**
 * Triangular face with three vertex indices. Vertex indices start with 0.
 */
struct Face {
    pub v1: usize,
    pub v2: usize,
    pub v3: usize,
}

struct Mesh {
    pub vertices: Vec<Vertex>,
    pub faces: Vec<Face>,
}

impl Mesh {
    fn cylinder() -> Self {
        let half_height = 5.0;
        let radius = 1.0;
        let linear_segments = 30;
        let radial_segments = 30;

        // Vertex indices: i * radial_segments + j
        let mut vertices: Vec<Vertex> = vec![];
        for i in 0..=linear_segments {
            let i_unipolar = (i as f32) / (linear_segments as f32);
            let i_bipolar = i_unipolar * 2.0 - 1.0;
            for j in 0..radial_segments {
                let theta = (j as f32) * std::f32::consts::TAU / (radial_segments as f32);
                let x = theta.cos() * radius;
                let y = theta.sin() * radius;
                let z = i_bipolar * half_height; 
                vertices.push(Vertex { x, y, z });
            }
        }

        let mut faces: Vec<Face> = vec![];
        for i in 0..linear_segments {
            for j in 0..radial_segments {
                let v1 = i * radial_segments + j;
                let v2 = i * radial_segments + (j + 1) % radial_segments;
                let v3 = (i + 1) * radial_segments + j;
                let v4 = (i + 1) * radial_segments + (j + 1) % radial_segments;

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

    fn transform_pipe(self, pipe: &Pipe, w: f32) -> Self {
        let new_vertices = self.vertices.into_iter().map(|vertex| {
            Vertex {
                x: vertex.x * pipe.a + vertex.y * pipe.b + pipe.c * w,
                y: -vertex.x * pipe.b + vertex.y * pipe.a + pipe.c * vertex.z,
                z: vertex.z,
            }
        }).collect::<Vec<_>>();
        Mesh { vertices: new_vertices, faces: self.faces }
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

    fn merge(meshes: Vec<Self>) -> Self {
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

    fn write_obj<W: Write>(&self, buffer: &mut W) -> std::io::Result<()> {
        for vertex in &self.vertices {
            write!(buffer, "v {} {} {}\n", vertex.x, vertex.y, vertex.z)?;
        }
        for face in &self.faces {
            // OBJ vertex indices start from 1.
            write!(buffer, "f {} {} {}\n", face.v1 + 1, face.v2 + 1, face.v3 + 1)?;
        }
        Ok(())
    }
}

fn main() -> std::io::Result<()> {
    let w = 0.0;
    let mesh = Mesh::merge(vec![
        Pipe { a: 1.0, b: 0.5, c: 0.2 },
        Pipe { a: 1.0, b: 0.0, c: 0.0 },
    ].iter().map(|pipe| {
        Mesh::cylinder().transform_pipe(&pipe, w)
    }).collect::<Vec<_>>());
    let mut buffer = File::create("out.obj")?;
    mesh.write_obj(&mut buffer)?;
    Ok(())
}
