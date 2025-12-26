#![allow(unused)]
extern crate approx;

use core::f64;
use std::fs;
use std::error::Error;
use std::path;
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
mod config;
use crate::config::Config;
use crate::polytwister::Polytwister;
use crate::utils::linspace;

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Args {
    /// Path to a JSON config file
    #[arg(short, long, value_name = "FILE")]
    config: Option<PathBuf>,
    
    #[command(subcommand)]
    command: Commands,
}


#[derive(Subcommand)]
enum Commands {
    Animation {
        /// Input JSON file describing the geometry of the polytwister.
        input_json: PathBuf,

        /// Number of animation frames.
        frames: usize,

        /// Output directory.
        output_dir: PathBuf,
    },
    Section {
        /// Input JSON file describing the geometry of the polytwister.
        input_json: PathBuf,

        /// W cross section
        w: f64,

        /// Output PLY mesh file.
        output_ply: PathBuf,
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();

    let config: Config = if let Some(config_path) = args.config {
        let mut string = String::new();
        let mut config_file = File::open(config_path)?;
        config_file.read_to_string(&mut string)?;
        serde_json::from_str(&string)?
    } else {
        Default::default()
    };

    match &args.command {
        Commands::Section { input_json, w, output_ply } => {
            let config: Config = Default::default();

            let mut string = String::new();
            let mut file = File::open(input_json)?;
            file.read_to_string(&mut string)?;
            let polytwister: Polytwister = serde_json::from_str(&string)?;

            let mesh = polytwister.as_colored_mesh(*w, &config);

            let mut buffer = File::create(output_ply)?;
            mesh.write_ply(&mut buffer)?;
            Ok(())
        },
        Commands::Animation { input_json, output_dir, frames } => {
            let config: Config = Default::default();

            let mut string = String::new();
            let mut file = File::open(input_json)?;
            file.read_to_string(&mut string)?;
            let polytwister: Polytwister = serde_json::from_str(&string)?;

            std::fs::create_dir(output_dir)?;

            for (i, w) in linspace(-1.0, 1.0, *frames).into_iter().enumerate() {
                eprintln!("frame {i}, w = {w}");
                let file_name = format!("section_{i:04}.ply");
                let output_ply = output_dir.join(file_name);
                let mesh = polytwister.as_colored_mesh(w, &config);
                let mut buffer = File::create(output_ply)?;
                mesh.write_ply(&mut buffer)?;
            }

            Ok(())
        }
    }
}