#![allow(unused)]
extern crate approx;

use core::f64;
use std::{any, fs};
use std::error::Error;
use std::path;
use std::path::PathBuf;
use std::{fs::File, io::Read};
use serde::Deserialize;
extern crate nalgebra as na;
use na::{Point3, Vector4};
use clap::{Parser, Subcommand, ValueEnum};

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
use crate::mesh::Mesh;

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
    /// Export a series of cross section meshes to a directory for an animation.
    /// 
    /// The meshes are exported in a Stanford PLY format, separated by frame and element type as
    /// follows:
    /// 
    /// ```
    /// out_dir/frame_0000_rings.ply
    /// out_dir/frame_0000_strips.ply
    /// out_dir/frame_0000_twisters_1.ply
    /// out_dir/frame_0000_twisters_2.ply
    /// out_dir/frame_0001_rings.ply
    /// out_dir/frame_0001_strips.ply
    /// out_dir/frame_0001_twisters_1.ply
    /// out_dir/frame_0001_twisters_2.ply
    /// ```
    /// 
    /// and so forth. twisters_1 and twisters_2 are the two orbits of the twisters. If the
    /// polytwister has only one orbit, all twisters_2 meshes will be empty.
    Animation {
        /// Input polytwister geometry file.
        /// 
        /// To get one of these files, use the "export-geometry" script in the Polytwisters JS app.
        input_json: PathBuf,

        /// Output directory.
        /// 
        /// The directory will be created. It is an error if the directory already exists.
        output_dir: PathBuf,

        /// Number of animation frames.
        /// 
        /// The W coordinates will be evenly spaced from -1 to +1.
        #[arg(short = 'n', long, default_value_t = 24)]
        frames: usize,
    },

    /// Export a single cross section of a polytwister as a mesh in the Stanford PLY format.
    /// 
    /// You can export the rings, strips, or individual twister orbits as meshes by providing the
    /// relevant options. You can also use the `--merged` option to generate a mesh that merges them
    /// all together with colors for visualization.
    Section {
        /// Input polytwister geometry file.
        /// 
        /// To get one of these files, use the "export-geometry" script in the Polytwisters JS app.
        input_json: PathBuf,

        /// W coordinate for the cross section.
        #[arg(short, default_value_t = 0.1)]
        w: f64,

        /// Output PLY mesh with everything: rings, strips, and twisters.
        /// 
        /// This is just for quickly loading the file to inspect in a mesh viewer, so there are not
        /// a lot of customization options here.
        #[arg(long = "merged")]
        merged_path: Option<PathBuf>,

        /// Output PLY mesh for ring cross sections.
        #[arg(long = "rings")]
        rings_path: Option<PathBuf>,

        /// Output PLY mesh for strip cross sections.
        #[arg(long = "strips")]
        strips_path: Option<PathBuf>,

        /// Output PLY mesh for cross sections of twisters in orbit 1.
        #[arg(long = "twisters-1")]
        twisters_path_1: Option<PathBuf>,

        /// Output PLY mesh for cross sections twisters in orbit 2.
        #[arg(long = "twisters-2")]
        twisters_path_2: Option<PathBuf>,
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
        Commands::Section { input_json, w, merged_path, rings_path, strips_path, twisters_path_1, twisters_path_2 } => {
            let config: Config = Default::default();

            let mut string = String::new();
            let mut file = File::open(input_json)?;
            file.read_to_string(&mut string)?;
            let polytwister: Polytwister = serde_json::from_str(&string)?;

            let mut any_output = false;
            let mesh = polytwister.as_meshes(*w, &config);
            if let Some(path) = merged_path {
                mesh.as_colored_mesh().write_ply_file(path);
                any_output = true;
            }
            if let Some(path) = rings_path {
                Mesh::merge(mesh.ring_meshes).write_ply_file(path);
                any_output = true;
            }
            if let Some(path) = strips_path {
                Mesh::merge(mesh.strip_meshes).write_ply_file(path);
                any_output = true;
            }
            if let Some(path) = twisters_path_1 {
                Mesh::merge(mesh.twister_meshes_orbit_1).write_ply_file(path);
                any_output = true;
            }
            if let Some(path) = twisters_path_2 {
                Mesh::merge(mesh.twister_meshes_orbit_2).write_ply_file(path);
                any_output = true;
            }
            if !any_output {
                eprintln!("Warning: no output mesh files provided. Try using --merged, --rings, --strips, --twisters-1, or --twisters-2.");
            }

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
                let prefix = format!("section_{i:04}");
                let mesh = polytwister.as_meshes(w, &config);
                mesh.write_plys(output_dir, &prefix);
            }

            Ok(())
        }
    }
}