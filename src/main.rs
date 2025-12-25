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
    let w = 0.1;
    let pipe_sections = vec![
        PipeSection::new(1.0, 0.0, 0.0, 0.0, w),
        PipeSection::new(0.0, 0.5, 0.5, 0.0, w),
        PipeSection::new(0.0, 0.0, 1.0, 0.0, w),
    ];

    let cylinder_options = CylinderMeshOptions {
        half_length: 5.0,
        linear_segments: 128,
        radial_segments: 128,
    };

    let mut meshes = vec![];

    for (i, pipe_section) in pipe_sections.iter().enumerate() {
        let mut mesh = pipe_section.as_mesh(&cylinder_options);
        for (j, pipe_section_2) in pipe_sections.iter().enumerate() {
            if i != j {
                mesh = mesh.filter_vertices(|p| pipe_section_2.contains(p));
            }
        }
        meshes.push(mesh);
    }

    let strip_thickness = 0.05;
    let strip_radial_resolution = 16;
    let strip_linear_resolution = 128;
    for (i, pipe_section_1) in pipe_sections.iter().enumerate() {
        for (j, pipe_section_2) in pipe_sections.iter().enumerate() {
            if i == j {
                continue;
            }

            let strip: Vec<Polyline> = pipe_section_1.as_cylinder().intersect_cylinder(&pipe_section_2.as_cylinder(), strip_linear_resolution);
            let strip_meshes = strip.iter().map(|polyline|
                polyline.as_mesh(strip_thickness, strip_radial_resolution)
            ).collect::<_>();
            let mut mesh = Mesh::merge(strip_meshes);
            for (k, pipe_section_3) in pipe_sections.iter().enumerate() {
                if k == i || k == j {
                    continue;
                }
                mesh = mesh.filter_vertices(|p| pipe_section_3.contains(p));
            }
            meshes.push(mesh);
        }
    }

    let mesh = Mesh::merge(meshes);

    let mut buffer = File::create("out.obj")?;
    mesh.write_obj(&mut buffer)?;
    Ok(())
}