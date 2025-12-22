#[macro_use]
extern crate approx;

use core::f64;
use std::{fs::File, io::Read};
use serde::Deserialize;
extern crate nalgebra as na;
use na::{Vector3};

mod strip_curves;
mod pipe_section;
mod mesh;
mod utils;
mod ellipse_spacing;
use crate::pipe_section::{PipeSection};
use crate::mesh::{Mesh};


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
        PipeSection { a: log[0], b: log[1], c: log[2], d: 0.0, w }
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