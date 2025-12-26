#![allow(unused)]
extern crate approx;

use core::f64;
use std::path::PathBuf;
use std::{fs::File, io::Read};
use serde::Deserialize;
extern crate nalgebra as na;
use na::{Point3, Vector4};
use clap::{Parser, Subcommand};

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
    config: Option<PathBuf>,
    
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    CrossSection {
        /// Input JSON file describing the geometry of the polytwister.
        input_json: PathBuf,

        /// W cross section
        w: f64,

        /// Output PLY mesh file.
        output_ply: PathBuf,
    }
}

fn main() -> std::io::Result<()> {
    let args = Args::parse();

    match &args.command {
        Commands::CrossSection { input_json, w, output_ply } => {
            let mut string = String::new();
            let mut file = File::open(input_json)?;
            file.read_to_string(&mut string)?;

            let polytwister: Polytwister = serde_json::from_str(&string)?;
            let mesh = polytwister.as_mesh(*w);

            let mut buffer = File::create(output_ply)?;
            mesh.write_ply(&mut buffer)?;
            Ok(())
        }
    }
}