use nalgebra as na;
use na::{Affine3, Point3};
use crate::mesh::{Mesh, Face};
use std::f64;

pub struct Curve {
    pub points: Vec<Point3<f64>>
}

impl Curve {
    fn point(&self, index: usize) -> Point3<f64> {
        self.points[index.rem_euclid(self.points.len())]
    }

    pub fn as_mesh(&self) -> Mesh {
        let num_points = self.points.len();
        let mut vertices = vec![];
        let radial_segments = 4;
        let thickness = 0.05;
        for (point_index, point) in self.points.iter().enumerate() {
            let prev = self.point(point_index - 1);
            let next = self.point(point_index + 1);
            let v_next = (next - point).normalize();
            let v_prev = (point - prev).normalize();
            let x = v_next.cross(&v_prev).normalize();
            let y = x.cross(&v_next).normalize();
            for radial_index in 0..radial_segments {
                let angle = radial_index as f64 / radial_segments as f64 * f64::consts::TAU;
                let cos = angle.cos();
                let sin = angle.sin();
                let mesh_point = point + (x * cos + y * sin) * thickness;
                vertices.push(mesh_point);
            }
        }

        let mut faces: Vec<Face> = vec![];

        for i in 0..num_points {
            for j in 0..radial_segments {
                let v1 = i * radial_segments + j;
                let v2 = i * radial_segments + (j + 1) % radial_segments;
                let v3 = (i + 1) % num_points * radial_segments + j;
                let v4 = (i + 1) % num_points * radial_segments + (j + 1) % radial_segments;

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

    pub fn transform(self, transform: &Affine3<f64>) -> Curve {
        let points = self.points.into_iter().map(|vertex|
            transform.transform_point(&vertex)
        ).collect::<Vec<_>>();
        Curve { points }
    }
}