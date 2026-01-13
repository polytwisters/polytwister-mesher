use nalgebra as na;
use na::{Affine3, Point3, Vector3};
use crate::mesh::{self, Face, Mesh, Vertex};
use core::num;
use std::f64;

/**
 * A polyline in 3D space.
 */
pub struct Polyline {
    pub points: Vec<Point3<f64>>
}

impl Polyline {
    fn plane(&self, i: usize) -> (Vector3<f64>, Vector3<f64>) {
        let i = i.clamp(1, self.points.len() - 2);
        let prev = self.points[i - 1];
        let point = self.points[i];
        let next = self.points[i + 1];
        let v_next = (next - point).normalize();
        let v_prev = (point - prev).normalize();
        let x = v_next.cross(&v_prev).normalize();
        let y = x.cross(&v_next).normalize();
        (x, y)
    }

    pub fn as_mesh(
        &self,
        thickness: f64,
        radial_segments: usize,
    ) -> Mesh {
        let num_points = self.points.len();

        // Doesn't make sense to do less than 3 points as cross products will be undefined.
        if num_points < 3 {
            return Mesh::empty();
        }

        let mut vertices = vec![];
        for i in 0..num_points {
            let point = self.points[i];
            let (x, y) = self.plane(i);
            for radial_index in 0..radial_segments {
                let angle = radial_index as f64 / radial_segments as f64 * f64::consts::TAU;
                let cos = angle.cos();
                let sin = angle.sin();
                let normal = x * cos + y * sin;
                let mesh_point = point + normal * thickness;
                let vertex = Vertex::new(mesh_point, normal);
                vertices.push(vertex);
            }
        }

        let mut faces: Vec<Face> = vec![];

        for i in 0..num_points - 1 {
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
                    v2: v3,
                    v3: v2,
                });
                faces.push(Face {
                    v1: v2,
                    v2: v3,
                    v3: v4,
                });
            }
        }

        Mesh { vertices, faces }
    }

    pub fn transform(self, transform: &Affine3<f64>) -> Polyline {
        let points = self.points.into_iter().map(|vertex|
            transform.transform_point(&vertex)
        ).collect::<Vec<_>>();
        Polyline { points }
    }
}