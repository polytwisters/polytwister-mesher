#![allow(unused)]
extern crate approx;

use core::f64;
use std::{fs::File, io::Read};
extern crate nalgebra as na;

mod cylinder_intersections;
mod pipe_section;
mod cylinder;
mod mesh;
mod utils;
mod ellipse_spacing;
mod polyline;
mod ring;
mod polytwister;
use crate::polytwister::Polytwister;

fn main() -> std::io::Result<()> {
    let mut string = String::new();
    let mut file = File::open("tetratwister.json")?;
    file.read_to_string(&mut string)?;

    let polytwister: Polytwister = serde_json::from_str(&string)?;
    let w = 0.1;
    let mesh = polytwister.as_mesh(w);

    let mut buffer = File::create("out.ply")?;
    mesh.write_ply(&mut buffer)?;
    Ok(())
}