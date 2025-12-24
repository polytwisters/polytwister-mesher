use nalgebra as na;
use na::{Affine3, Point3, Transform3};
use crate::mesh::Mesh;

pub struct Curve {
    pub points: Vec<Point3<f64>>
}

impl Curve {
    pub fn as_mesh(&self) -> Mesh {
        let meshes = self.points.iter().map(|point|
            Mesh::marker(&point)
        ).collect::<_>();
        Mesh::merge(meshes)
    }

    pub fn transform(self, transform: &Affine3<f64>) -> Curve {
        let points = self.points.into_iter().map(|vertex|
            transform.transform_point(&vertex)
        ).collect::<Vec<_>>();
        Curve { points }
    }
}