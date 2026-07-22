#![allow(unused)]
extern crate approx;

use core::f64;
use std::{any, fs};
use std::error::Error;
use std::path;
use std::path::PathBuf;
use std::{fs::File, io::Read, io::Write};
use serde::{Deserialize, Serialize};
extern crate nalgebra as na;
use na::{Point3, Vector4};
use clap::{Parser, Subcommand, ValueEnum};

use env_logger::Env;
use log::{warn, info};

mod cylinder_intersections;
mod pipe_section;
mod cylinder;
mod cylinder_curve;
mod mesh;
mod utils;
mod ellipse_spacing;
mod polyline;
mod ring;
mod uniform_polytwister;
mod c2;
mod convex_polytwister;
mod config;
mod marching_squares;
mod polytwister;
mod elements;

use crate::config::Config;
use crate::polytwister::Polytwister;
use crate::uniform_polytwister::{UniformPolytwister, PolytwisterDatabase};
use crate::convex_polytwister::{ConvexPolytwister, ConvexPolytwisterSpec};
use crate::utils::linspace;
use crate::mesh::{Mesh, MeshLike};


const MAX_FRAMES: usize = 10_000;
const DEFAULT_W: f64 = 0.1;

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Args {
    /// Path to a JSON config file
    #[arg(short, long, value_name = "FILE")]
    config: Option<PathBuf>,

    /// Name of the polytwister. You can use a full name, acronym, ID, or a JSON description of a
    /// convex polytwister.
    polytwister: String,

    /// Input polytwister database file.
    /// 
    /// A file is provided for you in the repo at "./polytwisters.json", which this defaults to,
    /// so you don't need to provide this option if your working directory contains that file.
    /// To get one of these files, use the "export-geometry" script in the Polytwisters JS app.
    #[arg(short = 'd', long)]
    database: Option<PathBuf>,

    /// W coordinate if exporting a single cross section.
    /// 
    /// It is an error to use -w and --frames together. If neither -w or --frames is specified, the
    /// default is equivalent to -w 0.1.
    #[arg(short)]
    w: Option<f64>,

    /// Number of animation frames. The W coordinates will be evenly spaced from -1 to +1 inclusive.
    ///
    /// It is an error to use -w and --frames together.
    #[arg(short = 'n', long)]
    frames: Option<usize>,

    /// Output a single merged mesh with all rings, strips, and twisters instead of a directory of
    /// meshes. This is used for quickly inspecting the result.
    /// 
    /// It is an error to use --merged and --frames together.
    #[arg(long)]
    merged: bool,

    /// Output path.
    output_path: PathBuf,
}

fn load_polytwister(polytwister_name: &String, database_path: &PathBuf) -> Result<Box<dyn Polytwister>, Box<dyn Error>> {
    if polytwister_name.starts_with("{") {
        let polytwister_spec: ConvexPolytwisterSpec = serde_json::from_str(polytwister_name)?;
        let polytwister = polytwister_spec.to_convex_polytwister().normalize();
        Ok(Box::new(polytwister))
    } else {
        let mut string = String::new();
        let mut file = File::open(database_path)?;
        file.read_to_string(&mut string)?;
        let polytwister_database: PolytwisterDatabase = serde_json::from_str(&string)?;
        let polytwister = polytwister_database.find(&polytwister_name)?.normalize();
        Ok(Box::new(polytwister))
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let env = Env::default()
        .filter_or("LOG_LEVEL", "info");
    env_logger::Builder::from_env(env).format_timestamp(None).init();

    let args = Args::parse();

    if let Some(_) = args.frames {
        if let Some(_) = args.w {
            return Err(String::from("-w and --frames cannot be used together").into());
        }
        if args.merged {
            return Err(String::from("--merged and --frames cannot be used together").into());
        }
    }

    let config: Config = if let Some(config_path) = args.config {
        let mut string = String::new();
        let mut config_file = File::open(config_path)?;
        config_file.read_to_string(&mut string)?;
        serde_json::from_str(&string)?
    } else {
        Default::default()
    };

    let polytwister_name = args.polytwister;
    let default_database_path = PathBuf::from("./polytwisters.json");
    let database_path = args.database.clone().unwrap_or(default_database_path);
    let polytwister = load_polytwister(&polytwister_name, &database_path)?;

    let out_path = args.output_path;

    if let Some(frames) = args.frames {
        if frames > MAX_FRAMES {
            return Err(String::from("Too many frames").into());
        }
        std::fs::create_dir(&out_path)?;
        let w_values = linspace(-1.0, 1.0, frames);
        for (i, w_ref) in w_values.iter().enumerate() {
            let w = *w_ref;
            info!("Frame {i}/{frames}, w = {w}");
            let section_dir = out_path.join(format!("section_{i:04}"));
            let mesh = polytwister.as_meshes(w, &config);
            std::fs::create_dir(&section_dir)?;
            mesh.write_plys(&section_dir);
        }

        let manifest_path = out_path.join("manifest.json");
        let manifest = AnimationManifest { w_values };
        let manifest_json = serde_json::to_string(&manifest)?;
        {
            let mut buffer = File::create(manifest_path)?;
            buffer.write(manifest_json.as_bytes())?;
        }
    } else {
        let w = args.w.unwrap_or(DEFAULT_W);
        let mesh = polytwister.as_meshes(w, &config);
        if args.merged {
            mesh.as_colored_mesh().write_ply_file_and_log(&out_path, "Merged colored mesh");
        } else {
            std::fs::create_dir(&out_path)?;
            mesh.write_plys(&out_path);
        }
    }
    Ok(())
}

#[derive(Serialize)]
struct AnimationManifest {
    w_values: Vec<f64>
}