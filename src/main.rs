extern crate approx;

use core::f64;
use std::{fs::File, io::Read};
use serde::Deserialize;
extern crate nalgebra as na;
use na::Point3;

mod cylinder_intersections;
mod pipe_section;
mod cylinder;
mod mesh;
mod utils;
mod ellipse_spacing;
mod polyline;
use crate::cylinder::{Cylinder, CylinderMeshOptions};
use crate::pipe_section::{PipeSection};
use crate::mesh::{Mesh};
use crate::polyline::Polyline;


#[derive(Deserialize)]
#[serde(rename_all="camelCase")]
struct Polytwister {
    logs: Vec<Vec<f64>>,
}

fn main() -> std::io::Result<()> {
    let cylinder = Cylinder::base();
    let cylinder_2 = Cylinder::example();

    let cylinder_options = CylinderMeshOptions {
        half_length: 5.0,
        linear_segments: 32,
        radial_segments: 32,
    };

    let cylinder_mesh = cylinder.as_mesh(&cylinder_options);
    let cylinder_mesh_2 = cylinder_2.as_mesh(&cylinder_options);

    let strip_thickness = 0.05;
    let strip_radial_resolution = 16;
    let strips: Vec<Polyline> = cylinder.intersect(&cylinder_2, 64);
    let strip_meshes = strips.iter().map(|polyline|
        polyline.as_mesh(strip_thickness, strip_radial_resolution)
    ).collect::<_>();
    let mesh = Mesh::merge(strip_meshes);
    let mesh = Mesh::merge(vec![cylinder_mesh, cylinder_mesh_2, mesh]);

    let mut buffer = File::create("out.obj")?;
    mesh.write_obj(&mut buffer)?;
    Ok(())
}