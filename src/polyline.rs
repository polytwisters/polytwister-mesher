use nalgebra as na;
use na::{Affine3, Point3};
use crate::mesh::{self, Face, Mesh, Vertex};
use core::num;
use std::f64;

/**
 * A polyline forming a closed loop in 3D space.
 */
pub struct Polyline {
    pub points: Vec<Point3<f64>>
}

impl Polyline {
    fn point(&self, index: isize) -> Point3<f64> {
        self.points[index.rem_euclid(self.points.len() as isize) as usize]
    }

    /**
     * Produce a mesh visualizing this polyline as a thick tube.
     */
    pub fn as_mesh(&self, thickness: f64, radial_segments: usize) -> Mesh {
        let num_points = self.points.len();

        // Doesn't make sense to do less than 3 points as cross products will be undefined.
        if num_points < 3 {
            return Mesh::empty();
        }

        let mut vertices = vec![];
        for (point_index, point) in self.points.iter().enumerate() {
            let prev = self.point(point_index as isize - 1);
            let next = self.point(point_index as isize + 1);
            let v_next = (next - point).normalize();
            let v_prev = (point - prev).normalize();
            let x = v_next.cross(&v_prev).normalize();
            let y = x.cross(&v_next).normalize();
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