use nalgebra as na;
use na::{Affine3, Point3, Vector3};
use crate::mesh::{self, Face, Mesh, Vertex};
use core::num;
use std::f64;

/**
 * A polyline in 3D space.
 */
pub struct Polyline {
    pub points: Vec<Point3<f64>>,
    pub closed: bool,
}

impl Polyline {
    fn num_points(&self) -> usize {
        self.points.len()
    }

    fn point(&self, i: isize) -> Point3<f64> {
        if self.closed {
            let index = i.rem_euclid(self.num_points() as isize) as usize;
            self.points[index]
        } else {
            self.points[i as usize]
        }
    }

    fn plane(&self, i: usize) -> (Vector3<f64>, Vector3<f64>) {
        // Index of the polygon we're using to estimate the curvature vector.
        let j = if self.closed {
            i
        } else {
            i.clamp(1, self.num_points() - 2)
        } as isize;
        let prev = self.point(j - 1);
        let point = self.point(j);
        let next = self.point(j + 1);
        let v_next = (next - point).normalize();
        let v_prev = (point - prev).normalize();
        let x = v_next.cross(&v_prev).normalize();
        let y = x.cross(&v_next).normalize();
        if !x.x.is_finite() {
            dbg!(prev, point, next, v_next, v_prev);
            panic!();
        }
        (x, y)
    }

    pub fn as_mesh(
        &self,
        thickness: f64,
        radial_segments: usize,
    ) -> Mesh {
        let num_points = self.num_points();

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

        let num_segments = if self.closed {
            num_points
        } else {
            // Leave off the last segment connecting the first and last points for an open polyline.
            num_points - 1
        };

        for i in 0..num_segments {
            for j in 0..radial_segments {
                let i1 = i;
                let i2 = (i + 1) % num_points;
                let j1 = j;
                let j2 = (j + 1) % radial_segments;
                let v1 = i1 * radial_segments + j1;
                let v2 = i1 * radial_segments + j2;
                let v3 = i2 * radial_segments + j1;
                let v4 = i2 * radial_segments + j2;

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
        let closed = self.closed;
        let points = self.points.into_iter().map(|vertex|
            transform.transform_point(&vertex)
        ).collect::<Vec<_>>();
        Polyline { points, closed }
    }
}