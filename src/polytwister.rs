use std::path::PathBuf;

use clap::Error;
use log;
use serde::Deserialize;
use na::{Point3, Vector4};
use crate::config::{Config, CylinderMeshConfig, RingMeshConfig, TorusMeshConfig};
use crate::cylinder::{Cylinder};
use crate::pipe_section::{PipeSection, StripSection, TwisterSection};
use crate::mesh::{Color, Mesh, ColoredMesh, MeshLike};
use crate::polyline::Polyline;
use crate::ring::{self, RingSection};

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
    geometry: Polytwister,
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
pub struct Polytwister {
    polyhedron: Polyhedron,
    pipes: Vec<Vector4<f64>>,
    orthogonal_pipes: Vec<Vector4<f64>>,
    rings: Vec<Vector4<f64>>,
    twister_fillings: Vec<Vec<FillingRegion>>,
    bloated: bool,
}

pub struct PolytwisterMeshes {
    pub ring_meshes: Vec<Mesh>,
    pub strip_meshes: Vec<Mesh>,
    pub twister_meshes_orbit_1: Vec<Mesh>,
    pub twister_meshes_orbit_2: Vec<Mesh>,
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

impl Polytwister {
    fn ring_cross_section(&self, index: usize, w: f64) -> RingSection {
        RingSection::from_vector4(&self.rings[index], w)
    }

    fn pipe_cross_section(&self, index: usize, w: f64) -> PipeSection {
        PipeSection::from_vector4(&self.pipes[index], w)
    }

    fn orthogonal_pipe_cross_section(&self, index: usize, w: f64) -> PipeSection {
        PipeSection::from_vector4(&self.orthogonal_pipes[index], w)
    }

    fn pipe_cross_sections(&self, w: f64) -> Vec<PipeSection> {
        (0..self.pipes.len()).map(|i| { self.pipe_cross_section(i, w) }).collect::<Vec<_>>()
    }

    fn orthogonal_pipe_cross_sections(&self, w: f64) -> Vec<PipeSection> {
        (0..self.pipes.len()).map(|i| { self.orthogonal_pipe_cross_section(i, w) }).collect::<Vec<_>>()
    }

    fn ring_cross_sections(&self, w: f64) -> Vec<RingSection> {
        (0..self.pipes.len()).map(|i| { self.ring_cross_section(i, w) }).collect::<Vec<_>>()
    }

    pub fn twister_sections(&self, w: f64, config: &Config) -> Vec<TwisterSection> {
        let pipe_sections = self.pipe_cross_sections(w);

        let mut twister_sections = vec![];
        for (pipe_index, pipe_section) in pipe_sections.iter().enumerate() {
            let face_orbit = self.polyhedron.faces[pipe_index].orbit;
            let filling = self.twister_fillings[face_orbit as usize].clone();
            let orthogonal_pipe_section = PipeSection::from_vector4(&self.orthogonal_pipes[pipe_index], w);
            let neighboring_pipe_sections = self.polyhedron.adjacent_face_indices(pipe_index)
                .iter().map(|pipe_index_2| pipe_sections[*pipe_index_2]).collect::<Vec<_>>();
            twister_sections.push(TwisterSection::new(
                pipe_section.clone(),
                orthogonal_pipe_section,
                neighboring_pipe_sections,
                filling
            ));
        }
        twister_sections
    }

    pub fn twister_orbit_as_meshes(&self, w: f64, orbit: u8, config: &Config) -> Vec<Mesh> {
        let twister_sections = self.twister_sections(w, config);
        let mut meshes = vec![];
        for (index, twister_section) in twister_sections.iter().enumerate() {
            if self.polyhedron.faces[index].orbit == orbit {
                let mesh = twister_section.as_mesh(&config.twisters);
                meshes.push(mesh);
            }
        }
        meshes
    }

    pub fn twister_orbit_as_mesh(&self, w: f64, orbit: u8, config: &Config) -> Mesh {
        Mesh::merge(self.twister_orbit_as_meshes(w, orbit, config))
    }

    pub fn rings_as_meshes(&self, w: f64, config: &Config) -> Vec<Mesh> {
        let ring_sections = self.ring_cross_sections(w);
        let mut meshes = vec![];
        for ring_section in ring_sections {
            let mesh = ring_section.as_mesh(&config.rings);
            meshes.push(mesh);
        }
        meshes
    }

    pub fn rings_as_mesh(&self, w: f64, config: &Config) -> Mesh {
        Mesh::merge(self.rings_as_meshes(w, config))
    }

    pub fn strips_as_meshes(&self, w: f64, config: &Config) -> Vec<Mesh> {
        let pipe_sections = self.pipe_cross_sections(w);
        let ring_sections = self.ring_cross_sections(w);
        let orthogonal_pipe_sections = self.orthogonal_pipe_cross_sections(w);

        self.polyhedron.edges.iter().enumerate().map(|(edge_index, edge)| {
            let adjacent_face_indices = self.polyhedron.edge_adjacent_face_indices(edge_index);
            if adjacent_face_indices.len() != 2 {
                panic!("Edge not adjacent to two faces");
            }
            let pipe_section_1 = pipe_sections[adjacent_face_indices[0]];
            let pipe_section_2 = pipe_sections[adjacent_face_indices[1]];
            let ring_section_1 = ring_sections[edge.vertex1];
            let ring_section_2 = ring_sections[edge.vertex2];
            let orthogonal_pipe_section = orthogonal_pipe_sections[adjacent_face_indices[0]];
            let torus_section = StripSection::new(
                &pipe_section_1,
                &pipe_section_2,
                &orthogonal_pipe_section,
                &ring_section_1,
                &ring_section_2,
                self.bloated
            );
            torus_section.as_mesh(&config.strips)
        }).collect::<Vec<_>>()
    }

    pub fn strips_as_mesh(&self, w: f64, config: &Config) -> Mesh {
        Mesh::merge(self.strips_as_meshes(w, config))
    }

    pub fn as_meshes(&self, w: f64, config: &Config) -> PolytwisterMeshes {
        PolytwisterMeshes {
            ring_meshes: self.rings_as_meshes(w, config),
            strip_meshes: self.strips_as_meshes(w, config),
            twister_meshes_orbit_1: self.twister_orbit_as_meshes(w, 0, config),
            twister_meshes_orbit_2:  self.twister_orbit_as_meshes(w, 1, config),
        }
    }

    pub fn as_colored_mesh(&self, w: f64, config: &Config) -> ColoredMesh {
        self.as_meshes(w, config).as_colored_mesh()
    }
}

impl PolytwisterMeshes {
    pub fn write_plys(&self, dir: &PathBuf, prefix: &str) -> std::io::Result<()> {
        let ring_mesh = Mesh::merge(self.ring_meshes.clone());
        let strip_mesh = Mesh::merge(self.strip_meshes.clone());
        let twister_mesh_orbit_1 = Mesh::merge(self.twister_meshes_orbit_1.clone());
        let twister_mesh_orbit_2 = Mesh::merge(self.twister_meshes_orbit_2.clone());

        ring_mesh.write_ply_file_and_log(&dir.join(format!("{prefix}_rings.ply")), "Ring mesh")?;
        strip_mesh.write_ply_file_and_log(&dir.join(format!("{prefix}_strips.ply")), "Strip mesh")?;
        twister_mesh_orbit_1.write_ply_file_and_log(&dir.join(format!("{prefix}_twisters_1.ply")), "Twister orbit 1 mesh")?;
        twister_mesh_orbit_2.write_ply_file_and_log(&dir.join(format!("{prefix}_twisters_2.ply")), "Twister orbit 2 mesh")?;

        Ok(())
    }

    /// Combine all ring, strip, and twister meshes into a single colored mesh.
    pub fn as_colored_mesh(&self) -> ColoredMesh {
        let strip_color = Color { red: 255, green: 255, blue: 255 };
        let ring_color = Color { red: 150, green: 150, blue: 150 };
        let twister_colors = vec![
            Color { red: 255, green: 0, blue: 238 },
            Color { red: 25, green: 25, blue: 255 },
        ];
        let ring_mesh = Mesh::merge(self.ring_meshes.clone());
        let strip_mesh = Mesh::merge(self.strip_meshes.clone());
        let twister_mesh_orbit_1 = Mesh::merge(self.twister_meshes_orbit_1.clone());
        let twister_mesh_orbit_2 = Mesh::merge(self.twister_meshes_orbit_2.clone());
        ColoredMesh {
            meshes: vec![
                (ring_mesh, ring_color),
                (strip_mesh, strip_color),
                (twister_mesh_orbit_1, twister_colors[0]),
                (twister_mesh_orbit_2, twister_colors[1]),
            ]
        }
    }
}

impl PolytwisterDatabase {
    /// Given a search string, find a matching polytwister. Matches by full name, acronym, ID,
    /// and index. Case-insensitive and _ and - may be substituted for spaces.
    pub fn find(&self, query: &str) -> Result<Polytwister, std::io::Error> {
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
            let pipe_section = polytwister.pipe_cross_section(face_index, w);
            let face = polytwister.polyhedron.faces[face_index].clone();
            for &vertex_index in face.vertices.iter() {
                let ring_section = polytwister.ring_cross_section(vertex_index, w);
                if let Some((point_1, point_2)) = ring_section.as_points() {
                    let tolerance = 1e-5;
                    assert!(pipe_section.boundary_contains(&point_1, tolerance));
                    assert!(pipe_section.boundary_contains(&point_2, tolerance));
                }
            }
        }
    }
}