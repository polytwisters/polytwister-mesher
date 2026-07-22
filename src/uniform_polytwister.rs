use std::io::pipe;
use std::path::PathBuf;

use clap::Error;
use log;
use serde::Deserialize;
use na::{Point3, Vector4};
use crate::config::{Config, CylinderMeshConfig, RingMeshConfig, TorusMeshConfig};
use crate::cylinder::{Cylinder};
use crate::pipe_section::{PipeSection, StripSection, TwisterSection};
use crate::mesh::{Color, Mesh, ColoredMesh, MeshLike, Polyline};
use crate::ring::{self, RingSection};
use crate::utils::torus_radius;
use crate::polytwister::{PolytwisterMeshes, Polytwister};

#[derive(Deserialize, Clone)]
#[serde(rename_all="camelCase")]
pub struct PolytwisterDatabase {
    polytwisters: Vec<PolytwisterWithDef>
}

#[derive(Deserialize, Clone)]
#[serde(rename_all="camelCase")]
pub struct PolytwisterDef {
    name: String,
    acronym: String,
    index: Option<u16>,
    symbol_string: String,
}

#[derive(Deserialize, Clone)]
#[serde(rename_all="camelCase")]
pub struct PolytwisterWithDef {
    geometry: UniformPolytwister,
    def: PolytwisterDef,
}

#[derive(Deserialize, Clone)]
#[serde(rename_all="camelCase")]
pub struct PolyhedronFace {
    vertices: Vec<usize>,
    edges: Vec<usize>,
    orbit: u8,
}

#[derive(Deserialize, Clone)]
#[serde(rename_all="camelCase")]
pub struct PolyhedronEdge {
    vertex1: usize,
    vertex2: usize,
}

#[derive(Deserialize, Clone)]
#[serde(rename_all="camelCase")]
pub struct Polyhedron {
    faces: Vec<PolyhedronFace>,
    edges: Vec<PolyhedronEdge>,
}

#[derive(Copy, Clone, Debug, Deserialize)]
pub enum RegionMode {
    Inner,
    Outer,
    Both
}

#[derive(Copy, Clone, Debug, Deserialize)]
#[serde(rename_all="camelCase")]
pub struct FillingRegion {
    pub order: u32,
    pub mode: RegionMode
}

#[derive(Deserialize, Clone)]
#[serde(rename_all="camelCase")]
pub struct UniformPolytwister {
    polyhedron: Polyhedron,
    pipes: Vec<Vector4<f64>>,
    orthogonal_pipes: Vec<Vector4<f64>>,
    rings: Vec<Vector4<f64>>,
    twister_fillings: Vec<Vec<FillingRegion>>,
    bloated: bool,
}

impl PolyhedronFace {
    fn adjacent_to_face(&self, other: &PolyhedronFace) -> bool {
        self.edges.iter().any(|e| other.edges.contains(e))
    }
}

impl Polyhedron {
    fn adjacent_face_indices(&self, face_index: usize) -> Vec<usize> {
        self.faces.iter().enumerate().filter_map(|(i, face)| {
            if i != face_index && face.adjacent_to_face(&self.faces[face_index]) {
                Some(i)
            } else {
                None
            }
        }).collect::<_>()
    }

    fn edge_adjacent_face_indices(&self, edge_index: usize) -> Vec<usize> {
        self.faces.iter().enumerate().filter_map(|(i, face)| {
            if face.edges.contains(&edge_index) {
                Some(i)
            } else {
                None
            }
        }).collect::<_>()
    }
}

impl UniformPolytwister {
    /// Return a new polytwister which is geometrically scaled by the given ratio.
    fn scale(&self, ratio: f64) -> Self {
        return UniformPolytwister {
            polyhedron: self.polyhedron.clone(),
            pipes: self.pipes.iter().map(|vec| vec / ratio).collect::<_>(),
            orthogonal_pipes: self.orthogonal_pipes.iter().map(|vec| vec / ratio).collect::<_>(),
            rings: self.rings.iter().map(|vec| vec * ratio).collect::<_>(),
            twister_fillings: self.twister_fillings.clone(),
            bloated: self.bloated,
        }
    }

    /// Return the polytwister's maximum distance from the origin.
    fn radius(&self) -> f64 {
        if self.bloated {
            let index = 0;
            let index_2 = self.polyhedron.adjacent_face_indices(index)[0];
            let pipe_1 = self.pipes[index];
            let pipe_2 = self.pipes[index_2];
            torus_radius(&pipe_1, &pipe_2)
        } else {
            self.rings[0].norm()
        }
    }

    /// Return a scaled version of this polytwister so that its radius is 1.0.
    pub fn normalize(&self) -> Self {
        let scale = 1.0 / self.radius();
        self.scale(scale)
    }

    fn ring_section(&self, index: usize, w: f64) -> RingSection {
        RingSection::from_vector4(&self.rings[index], w)
    }

    fn pipe_section(&self, index: usize, w: f64) -> PipeSection {
        PipeSection::from_vector4(&self.pipes[index], w)
    }

    fn orthogonal_pipe_section(&self, index: usize, w: f64) -> PipeSection {
        PipeSection::from_vector4(&self.orthogonal_pipes[index], w)
    }

    pub fn twister_section(&self, index: usize, w: f64) -> TwisterSection {
        let face_orbit = self.polyhedron.faces[index].orbit;
        let filling = self.twister_fillings[face_orbit as usize].clone();
        let pipe_section = self.pipe_section(index, w);
        let orthogonal_pipe_section = self.orthogonal_pipe_section(index, w);
        let neighboring_pipe_sections = self.polyhedron.adjacent_face_indices(index)
            .iter().map(|index_2| self.pipe_section(*index_2, w)
        ).collect::<Vec<_>>();
        TwisterSection::new(
            pipe_section,
            orthogonal_pipe_section,
            neighboring_pipe_sections,
            filling
        )
    }

    fn pipe_sections(&self, w: f64) -> Vec<PipeSection> {
        (0..self.pipes.len()).map(|i| { self.pipe_section(i, w) }).collect::<Vec<_>>()
    }

    fn orthogonal_pipe_sections(&self, w: f64) -> Vec<PipeSection> {
        (0..self.pipes.len()).map(|i| { self.orthogonal_pipe_section(i, w) }).collect::<Vec<_>>()
    }

    fn ring_sections(&self, w: f64) -> Vec<RingSection> {
        (0..self.rings.len()).map(|i| { self.ring_section(i, w) }).collect::<Vec<_>>()
    }

    pub fn twister_sections(&self, w: f64) -> Vec<TwisterSection> {
        (0..self.pipes.len()).map(|i| { self.twister_section(i, w) }).collect::<Vec<_>>()
    }

    pub fn twister_orbit_as_mesh(&self, w: f64, orbit: u8, config: &Config) -> Mesh {
        Mesh::merge(self.twister_orbit_as_meshes(w, orbit, config))
    }

    pub fn rings_as_mesh(&self, w: f64, config: &Config) -> Mesh {
        Mesh::merge(self.rings_as_meshes(w, config))
    }

    pub fn strip_section(&self, index: usize, w: f64) -> StripSection {
        let edge = &self.polyhedron.edges[index];
        let adjacent_face_indices = self.polyhedron.edge_adjacent_face_indices(index);
        if adjacent_face_indices.len() != 2 {
            panic!("Edge not adjacent to two faces");
        }
        let pipe_section_1 = self.pipe_section(adjacent_face_indices[0], w);
        let pipe_section_2 = self.pipe_section(adjacent_face_indices[1], w);
        let ring_section_1 = self.ring_section(edge.vertex1, w);
        let ring_section_2 = self.ring_section(edge.vertex2, w);
        let orthogonal_pipe_section = self.orthogonal_pipe_section(adjacent_face_indices[0], w);
        StripSection::new(
            &pipe_section_1,
            &pipe_section_2,
            &orthogonal_pipe_section,
            &ring_section_1,
            &ring_section_2,
            self.bloated
        )
    }

    pub fn strips_as_mesh(&self, w: f64, config: &Config) -> Mesh {
        Mesh::merge(self.strips_as_meshes(w, config))
    }
}

impl Polytwister for UniformPolytwister {
    fn twister_orbit_as_meshes(&self, w: f64, orbit: u8, config: &Config) -> Vec<Mesh> {
        let twister_sections = self.twister_sections(w);
        let mut meshes = vec![];
        for (index, twister_section) in twister_sections.iter().enumerate() {
            if self.polyhedron.faces[index].orbit == orbit {
                let mesh = twister_section.as_mesh(&config.twisters);
                meshes.push(mesh);
            }
        }
        meshes
    }

    fn rings_as_meshes(&self, w: f64, config: &Config) -> Vec<Mesh> {
        let ring_sections = self.ring_sections(w);
        let mut meshes = vec![];
        for ring_section in ring_sections {
            let mesh = ring_section.as_mesh(&config.rings);
            meshes.push(mesh);
        }
        meshes
    }

    fn strips_as_meshes(&self, w: f64, config: &Config) -> Vec<Mesh> {
        (0..self.polyhedron.edges.len()).map(|edge_index| {
            self.strip_section(edge_index, w).as_mesh(&config.strips)
        }).collect::<Vec<_>>()
    }
}

impl PolytwisterDatabase {
    /// Given a search string, find a matching polytwister. Matches by full name, acronym, ID,
    /// and index. Case-insensitive and _ and - may be substituted for spaces.
    pub fn find(&self, query: &str) -> Result<UniformPolytwister, std::io::Error> {
        let query = query.to_ascii_lowercase().replace("_", " ").replace("-", " ");
        for polytwister_with_def in self.polytwisters.iter() {
            let def = &polytwister_with_def.def;
            if (
                def.acronym == query
                || def.name == query
                || def.symbol_string == query
                || if let Some(index) = def.index { index.to_string() == query } else { false }
            ) {
                return Ok(polytwister_with_def.geometry.clone());
            }
        }
        Err(std::io::Error::other("Couldn't find polytwister."))
    }
}

#[cfg(test)]
mod test {
    use crate::ring::RingSectionResult;

    use super::*;
    use std::fs::File;
    use std::io::Read;

    /// Take a cross section of the cube twister and check that all ring cross sections are on their
    /// incident pipe cross sections.
    #[test]
    fn test_integration() {
        let mut string = String::new();
        let mut file = File::open("polytwisters.json").unwrap();
        file.read_to_string(&mut string).unwrap();
        let polytwister_database: PolytwisterDatabase = serde_json::from_str(&string).unwrap();
        let polytwister = polytwister_database.find("cubiter").unwrap();

        let w = 0.1;
        for (face_index, face) in polytwister.polyhedron.faces.iter().enumerate() {
            let pipe_section = polytwister.pipe_section(face_index, w);
            let face = polytwister.polyhedron.faces[face_index].clone();
            for &vertex_index in face.vertices.iter() {
                let ring_section = polytwister.ring_section(vertex_index, w);
                if let RingSectionResult::Points(point_1, point_2) = ring_section.as_points() {
                    let tolerance = 1e-5;
                    assert!(pipe_section.boundary_contains(&point_1, tolerance));
                    assert!(pipe_section.boundary_contains(&point_2, tolerance));
                }
            }
        }
    }
}