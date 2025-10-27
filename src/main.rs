#[macro_use]
extern crate approx;

use std::{fs::File, io::{Read, Write}};
use serde::Deserialize;
extern crate nalgebra as na;
use na::{Vector3, Vector2, Matrix3, Rotation3};

fn squared(x: f64) -> f64 {
    x * x
}

/**
 * A 3D cross section of a pipe. (PipeCrossSection felt too long.)
 */
#[derive(Clone, Copy, Debug)]
pub struct PipeSection {
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub w: f64,
}

/**
 * Triangular face with three vertex indices. Vector3<f64> indices start with 0.
 */
#[derive(Clone, Copy, Debug)]
struct Face {
    pub v1: usize,
    pub v2: usize,
    pub v3: usize,
}

struct Mesh {
    pub vertices: Vec<Vector3<f64>>,
    pub faces: Vec<Face>,
}

impl PipeSection {
    fn as_mesh(&self) -> Mesh {
        if self.a == 0.0 && self.b == 0.0 {
            let tmp = 1.0 / squared(self.c) - squared(self.w);
            if tmp <= 0.0 {
                return Mesh::empty();
            }
            let z = tmp.sqrt();
            return Mesh::merge(vec![
                Mesh::plane(z),
                Mesh::plane(-z),
            ])
        }
        Mesh::cylinder().transform_pipe(&self)
    }

    /**
     * Evaluate the scalar field associated with this pipe cross section:
     * 
     *     F(x, y, z) = (ax + by + cz)^2 + (bx - ay - cw)^2 - 1
     * 
     * The pipe cross section is given by the isosurface F(x, y, z) = 0. For other points, if
     * F(x, y, z) < 0 then the point is inside the pipe, F(x, y, z) > 0 is outside the pipe, and
     * F(x, y, z) = -1 is on the pipe's symmetry axis (plane of symmetry if a = b = 0).
     */
    fn scalar_field(&self, point: &Vector3<f64>) -> f64 {
        squared(self.a * point.x + self.b * point.y + self.c * point.z)
        + squared(self.b * point.x - self.a * point.y - self.c * self.w)
        - 1.0
    }

    fn abcw(&self) -> (f64, f64, f64, f64) {
        (self.a, self.b, self.c, self.w)
    }

    /**
     * Return an explicit parametrization of the line which is this pipe section's symmetry axis.
     * The parametrization is v(t) = v_0 + d * t and returned as (v_0, d) so that v_0 is the
     * starting point and d is the direction vector. d is always a unit vector.
     */
    fn axis_line(&self) -> (Vector3<f64>, Vector3<f64>) {
        let (a, b, c, w) = self.abcw();
        let tmp = 1.0 / (squared(a) + squared(b));
        let direction_vector = Vector3::new( -a * c * tmp, -b * c * tmp, 1.0).normalize();
        let start = Vector3::new(b * c * w * tmp, -a * c * w * tmp, 0.0);
        return (start, direction_vector);
    }

    /**
     * Return a 3x3 matrix that turns this pipe into the "base cylinder" x^2 + y^2 = 1, z in R.
     * It is assumed that a = b = 0 does not hold. Translation due to w is ignored.
     * 
     * This is called the "basic" backward matrix because it comes directly from the implicit
     * equations. Although it transforms the z-axis to the cylinder's axis of symmetry, it does not
     * preserve distance along that line.
     */
    fn basic_backward_matrix(&self) -> Matrix3<f64> {
        let (a, b, c, _) = self.abcw();
        Matrix3::new(
            a, b, c,
            b, -a, 0.0,
            0.0, 0.0, 1.0,
        )
    }

    /**
     * Return a 3x3 matrix that turns the "base cylinder" x^2 + y^2 = 1, z in R into this pipe.
     * It is assumed that a = b = 0 does not hold. Translation due to w is ignored.
     * 
     * This is the inverse of the basic_backward_matrix.
     */
    fn basic_forward_matrix(&self) -> Matrix3<f64> {
        let (a, b, c, _) = self.abcw();
        let tmp = 1.0 / (squared(a) + squared(b));
        Matrix3::new(
            a * tmp, b * tmp, -a * c * tmp,
            b * tmp, -a * tmp, -b * c * tmp,
            0.0, 0.0, 1.0,
        )
    }

    /**
     * Compute the two displacement vectors that span the elliptic cross-section of the pipe.
     */
    fn ellipse_vertex_displacements(&self) -> (Vector3<f64>, Vector3<f64>) {
        let (start, direction) = self.axis_line();
        let z = Vector3::z();
        // Matrix "r", when applied to the pipe cross section, has the effect of uprighting it so
        // that its symmetry axis is the z-axis.
        let r = Rotation3::rotation_between(&direction, &z).unwrap();
        let r_inv = Rotation3::rotation_between(&z, &direction).unwrap();
        let forward_ellipse: Matrix3<f64> = r * self.basic_forward_matrix();
        let displacement_1 = forward_ellipse * Vector3::x();
        let displacement_2 = forward_ellipse * Vector3::y();
        let point_a = r_inv * displacement_1;
        let point_b = r_inv * displacement_2;
        (start + point_a, start + point_b)
    }
}

impl Mesh {
    fn empty() -> Self {
        Mesh { vertices: vec![], faces: vec![] }
    }

    fn cylinder() -> Self {
        let half_height = 2.0;
        let radius = 1.0;
        let linear_segments = 500;
        let radial_segments = 200;

        // Vertex indices: i * radial_segments + j
        let mut vertices: Vec<Vector3<f64>> = vec![];
        for i in 0..=linear_segments {
            let i_unipolar = (i as f64) / (linear_segments as f64);
            let i_bipolar = i_unipolar * 2.0 - 1.0;
            for j in 0..radial_segments {
                let theta = (j as f64) * std::f64::consts::TAU / (radial_segments as f64);
                let x = theta.cos() * radius;
                let y = theta.sin() * radius;
                let z = i_bipolar * half_height; 
                vertices.push(Vector3::new(x, y, z));
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

    /**
     * Make a plane parallel to the xy-plane at coordinate z.
     */
    fn plane(z: f64) -> Self {
        let radius = 5.0;
        let segments = 30;

        // Vector3<f64> indices: i * segments + j
        let mut vertices: Vec<Vector3<f64>> = vec![];
        for i in 0..=segments {
            let i_unipolar = (i as f64) / (segments as f64);
            let i_bipolar = i_unipolar * 2.0 - 1.0;
            for j in 0..=segments {
                let j_unipolar = (j as f64) / (segments as f64);
                let j_bipolar = j_unipolar * 2.0 - 1.0;
                let x = i_bipolar * radius;
                let y = j_bipolar * radius;
                vertices.push(Vector3::new(x, y, z));
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

    fn transform_pipe(self, pipe: &PipeSection) -> Self {
        let new_vertices = self.vertices.into_iter().map(|vertex| {
            let (xp, yp, zp) = (vertex.x, vertex.y, vertex.z);
            let (a, b, c, w) = (pipe.a, pipe.b, pipe.c, pipe.w);
            /*
            Derivation: the pipe is P(a + bi, c + 0i) and we take its cross section at w. From the
            formula for pipes we have the implicit equation

                (ax + by + cz)^2 + (bx - ay - cw)^2 = 1.

            Let
                
                x' = ax + by + cz    (1)
                y' = bx - ay - cw
                z' = z
                
            so that x'^2 + y'^2 = 1 forms the implicit equation for the cylinder (x', y', z') of
            radius 1 and parallel to the z-axis. Thus we have an affine transformation relating
            (x, y, z) to (x', y', z').  Inverting this transformation lets us transform a base
            cylinder to a pipe cross section. The inversion is done by solving the system of
            equations (1) for x and y.
            */
            let tmp = 1.0 / (a * a + b * b);
            let z = zp;
            let x = (a * xp + b * yp + b * c * w - a * c * z) * tmp;
            let y = (b * xp - a * yp + a * c * w - b * c * z) * tmp;
            Vector3::new(x, y, z)
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

    /**
     * Given a predicate on vertex locations, return a new Mesh that removes all vertices that do
     * not satisfy that predicate, and any faces that are connected to said vertices.
     */
    fn filter_vertices<F: Fn(&Vector3<f64>) -> bool>(&self, predicate: F) -> Self {
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

#[derive(Deserialize)]
#[serde(rename_all="camelCase")]
struct Polytwister {
    logs: Vec<Vec<f64>>,
}

fn main() -> std::io::Result<()> {
    let mut string = String::new();
    let mut file = File::open("quasitetratwister.json")?;
    file.read_to_string(&mut string)?;

    let result: Polytwister = serde_json::from_str(&string)?;

    let w = 0.0;
    let pipes: Vec<PipeSection> = result.logs.iter().map(|log: &Vec<f64>| {
        PipeSection { a: log[0], b: log[1], c: log[2], w }
    }).collect::<Vec<_>>();

    let meshes = pipes.iter().enumerate().map(|(i, pipe)| {
        pipe.as_mesh().filter_vertices(|vertex: &Vector3<f64>| -> bool {
            for (j, pipe2) in pipes.iter().enumerate() {
                if i != j && (pipe2.scalar_field(vertex) >= 0.0) {
                    return false;
                }
            }
            true
        })
    }).collect::<Vec<_>>();

    let mesh = Mesh::merge(meshes);
    let mut buffer = File::create("out.obj")?;
    mesh.write_obj(&mut buffer)?;
    Ok(())
}

#[cfg(test)]
mod test {
    use na::Matrix3;

    use crate::PipeSection;

    #[test]
    fn test_pipe_axis() {
        let pipe = PipeSection { a: 1.2, b: -0.5, c: 0.4, w: 0.1 };
        let (start, direction) = pipe.axis_line();
        assert_abs_diff_eq!(pipe.scalar_field(&start), -1.0);
        let point_2 = start + direction;
        assert_abs_diff_eq!(pipe.scalar_field(&point_2), -1.0);
        let point_3 = start + direction * 1.2345;
        assert_abs_diff_eq!(pipe.scalar_field(&point_3), -1.0);
    }

    #[test]
    fn test_forward_inverse_matrices() {
        let pipe = PipeSection { a: 1.2, b: -0.5, c: 0.4, w: 0.1 };
        assert_abs_diff_eq!(
            pipe.basic_forward_matrix() * pipe.basic_backward_matrix(),
            Matrix3::identity()
        )
    }

    #[test]
    fn test_ellipse() {
        let pipe = PipeSection { a: 1.2, b: -0.5, c: 0.4, w: 0.1 };
        let (point_a, point_b) = pipe.ellipse_vertex_displacements();
        assert_abs_diff_eq!(pipe.scalar_field(&point_a), 0.0);
        assert_abs_diff_eq!(pipe.scalar_field(&point_b), 0.0);
    }
}