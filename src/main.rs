#![allow(unused)]
extern crate approx;

use core::f64;
use std::{fs::File, io::Read};
use serde::Deserialize;
extern crate nalgebra as na;
use na::{Point3, Vector4};
use serde::de::Visitor;

mod cylinder_intersections;
mod pipe_section;
mod cylinder;
mod mesh;
mod utils;
mod ellipse_spacing;
mod polyline;
mod ring;
use crate::cylinder::{Cylinder, CylinderMeshOptions};
use crate::pipe_section::{PipeSection, TorusMeshOptions, TorusSection, TwisterSection};
use crate::mesh::{Color, Mesh, MeshCollection};
use crate::polyline::Polyline;
use crate::ring::{RingSection, RingMeshOptions};

#[derive(Deserialize)]
#[serde(rename_all="camelCase")]
struct PolyhedronFace {
    vertices: Vec<usize>,
    edges: Vec<usize>,
    orbit: u8,
}

impl PolyhedronFace {
    fn adjacent_to(&self, other: &PolyhedronFace) -> bool {
        self.edges.iter().any(|e| other.edges.contains(e))
    }
}

#[derive(Deserialize)]
#[serde(rename_all="camelCase")]
struct Polyhedron {
    faces: Vec<PolyhedronFace>,
}

impl Polyhedron {
    fn adjacent_face_indices(&self, face_index: usize) -> Vec<usize> {
        self.faces.iter().enumerate().filter_map(|(i, face)| {
            if i != face_index && face.adjacent_to(&self.faces[face_index]) {
                Some(i)
            } else {
                None
            }
        }).collect::<_>()
    }
}


#[derive(Deserialize)]
#[serde(rename_all="camelCase")]
struct Polytwister {
    polyhedron: Polyhedron,
    logs: Vec<Vector4<f64>>,
    rings: Vec<Vector4<f64>>,
}

impl Polytwister {
    fn pipe_cross_sections(&self, w: f64) -> Vec<PipeSection> {
        self.logs.iter().map(|log: &Vector4<f64>| {
            PipeSection { a: log[0], b: log[1], c: log[2], d: 0.0, w }
        }).collect::<Vec<_>>()
    }

    fn ring_cross_sections(&self, w: f64) -> Vec<RingSection> {
        self.rings.iter().map(|ring: &Vector4<f64>| {
            RingSection { a: ring[0], b: ring[1], c: ring[2], d: 0.0, w }
        }).collect::<Vec<_>>()
    }

    fn as_mesh(&self, w: f64) -> MeshCollection {
        let pipe_sections = self.pipe_cross_sections(w);
        let ring_sections = self.ring_cross_sections(w);

        let cylinder_options = CylinderMeshOptions {
            half_length: 5.0,
            linear_segments: 128 * 2,
            radial_segments: 128 * 2,
        };
        let torus_options = TorusMeshOptions {
            thickness: 0.05,
            radial_segments: 16,
            linear_segments: 128,
        };
        let ring_options = RingMeshOptions {
            radius: 0.08,
            segments: 16,
            rings: 32, 
        };
        let twister_color = Color { red: 255, green: 0, blue: 0 };
        let strip_color = Color { red: 255, green: 255, blue: 255 };
        let ring_color = Color { red: 150, green: 150, blue: 150 };

        let mut twister_sections = vec![];
        for (pipe_index, pipe_section) in pipe_sections.iter().enumerate() {
            let neighboring_pipe_sections = self.polyhedron.adjacent_face_indices(pipe_index)
                .iter().map(|pipe_index_2| pipe_sections[*pipe_index_2]).collect::<Vec<_>>();
            twister_sections.push(TwisterSection {
                pipe_section: pipe_section.clone(),
                neighboring_pipe_sections
            });
        }

        let mut meshes = vec![];

        for twister_section in twister_sections {
            let mesh = twister_section.as_mesh(&cylinder_options);
            meshes.push((mesh, twister_color));
        }

        // Produce strip sections.
        for (i, pipe_section_1) in pipe_sections.iter().enumerate() {
            for (j, pipe_section_2) in pipe_sections.iter().enumerate() {
                // i == j is intersecting a pipe section with itself.
                // Ignoring i < j prevents doubling up strips since intersection is commutative.
                if i <= j {
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
                meshes.push((mesh, strip_color));
            }
        }

        for ring_section in ring_sections {
            let mesh = ring_section.as_mesh(&ring_options);
            meshes.push((mesh, ring_color));
        }

        MeshCollection { meshes }
    }
}

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