use nalgebra as na;
use na::Point3;
use crate::mesh::Mesh;

pub struct Curve {
    pub points: Vec<Point3<f64>>
}

impl Curve {

    fn as_mesh(&self) -> Mesh {
        Mesh::empty()
    }

}