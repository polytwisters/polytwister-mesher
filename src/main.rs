#![allow(unused)]
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
use crate::pipe_section::{PipeSection, TorusSection, TorusMeshOptions};
use crate::mesh::{Mesh};
use crate::polyline::Polyline;


#[derive(Deserialize)]
#[serde(rename_all="camelCase")]
struct PolytwisterJSON {
    logs: Vec<Vec<f64>>,
}

fn make_polytwister_mesh(pipe_sections: &Vec<PipeSection>) -> Mesh {
    let cylinder_options = CylinderMeshOptions {
        half_length: 5.0,
        linear_segments: 128,
        radial_segments: 128,
    };
    let torus_options = TorusMeshOptions {
        thickness: 0.05,
        radial_segments: 16,
        linear_segments: 128,
    };

    let mut meshes = vec![];

    // Produce twister sections.
    for (i, pipe_section) in pipe_sections.iter().enumerate() {
        let mut mesh = pipe_section.as_mesh(&cylinder_options);
        for (j, pipe_section_2) in pipe_sections.iter().enumerate() {
            if i != j {
                mesh = mesh.filter_vertices(|p| pipe_section_2.contains(p));
            }
        }
        meshes.push(mesh);
    }

    // Produce strip sections.
    for (i, pipe_section_1) in pipe_sections.iter().enumerate() {
        for (j, pipe_section_2) in pipe_sections.iter().enumerate() {
            if i == j {
                continue;
            }
            let torus_section = TorusSection::new(pipe_section_1, pipe_section_2);
            let mut mesh = torus_section.as_mesh(&torus_options);
            for (k, pipe_section_3) in pipe_sections.iter().enumerate() {
                if k == i || k == j {
                    continue;
                }
                mesh = mesh.filter_vertices(|p| pipe_section_3.contains(p));
            }
            meshes.push(mesh);
        }
    }

    Mesh::merge(meshes)
}

fn main() -> std::io::Result<()> {
    let mut string = String::new();
    let mut file = File::open("tetratwister.json")?;
    file.read_to_string(&mut string)?;

    let result: PolytwisterJSON = serde_json::from_str(&string)?;
    let w = 0.1;

    let pipe_sections: Vec<PipeSection> = result.logs.iter().map(|log: &Vec<f64>| {
        PipeSection { a: log[0], b: log[1], c: log[2], d: 0.0, w }
    }).collect::<Vec<_>>();

    let mesh = make_polytwister_mesh(&pipe_sections);

    let mut buffer = File::create("out.obj")?;
    mesh.write_obj(&mut buffer)?;
    Ok(())
}