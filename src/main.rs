#![allow(unused)]
extern crate approx;

use core::f64;
use std::path::PathBuf;
use std::{fs::File, io::Read};
use serde::Deserialize;
extern crate nalgebra as na;
use na::{Point3, Vector4};
use clap::Parser;

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


#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Args {
    /// Input JSON file describing the geometry of the polytwister.
    input_json: PathBuf,

    /// Output PLY mesh file.
    output_ply: PathBuf,

    config: Option<PathBuf>,
}

fn main() -> std::io::Result<()> {
    let args = Args::parse();

    let mut string = String::new();
    let mut file = File::open(args.input_json)?;
    file.read_to_string(&mut string)?;

    let polytwister: Polytwister = serde_json::from_str(&string)?;
    let w = 0.1;
    let mesh = polytwister.as_mesh(w);

    let mut buffer = File::create(args.output_ply)?;
    mesh.write_ply(&mut buffer)?;
    Ok(())
}