use serde::Deserialize;
use na::{Point3, Vector4};
use crate::config::{Config, CylinderMeshConfig, RingMeshConfig, TorusMeshConfig};
use crate::cylinder::{Cylinder};
use crate::pipe_section::{PipeSection, StripSection, TwisterSection};
use crate::mesh::{Color, Mesh, ColoredMeshCollection};
use crate::polyline::Polyline;
use crate::ring::{RingSection};

#[derive(Deserialize)]
#[serde(rename_all="camelCase")]
pub struct PolyhedronFace {
    vertices: Vec<usize>,
    edges: Vec<usize>,
    orbit: u8,
}

#[derive(Deserialize)]
#[serde(rename_all="camelCase")]
pub struct PolyhedronEdge {
    vertex1: usize,
    vertex2: usize,
}

#[derive(Deserialize)]
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

#[derive(Deserialize)]
#[serde(rename_all="camelCase")]
pub struct Polytwister {
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

impl Polytwister {
    fn pipe_cross_sections(&self, w: f64) -> Vec<PipeSection> {
        self.pipes.iter().map(|pipe: &Vector4<f64>| {
            PipeSection::from_vector4(&pipe, w)
        }).collect::<Vec<_>>()
    }

    fn orthogonal_pipe_cross_sections(&self, w: f64) -> Vec<PipeSection> {
        self.orthogonal_pipes.iter().map(|pipe: &Vector4<f64>| {
            PipeSection::from_vector4(&pipe, w)
        }).collect::<Vec<_>>()
    }

    fn ring_cross_sections(&self, w: f64) -> Vec<RingSection> {
        self.rings.iter().map(|ring: &Vector4<f64>| {
            RingSection::from_vector4(&ring, w)
        }).collect::<Vec<_>>()
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
            twister_sections.push(TwisterSection {
                pipe_section: pipe_section.clone(),
                orthogonal_pipe_section,
                neighboring_pipe_sections,
                filling
            });
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
        let orthogonal_pipe_sections = self.orthogonal_pipe_cross_sections(w);

        self.polyhedron.edges.iter().enumerate().map(|(edge_index, edge)| {
            let adjacent_face_indices = self.polyhedron.edge_adjacent_face_indices(edge_index);
            if adjacent_face_indices.len() != 2 {
                panic!("Edge not adjacent to two faces");
            }
            let pipe_section_1 = pipe_sections[adjacent_face_indices[0]];
            let pipe_section_2 = pipe_sections[adjacent_face_indices[1]];
            let orthogonal_pipe_section = orthogonal_pipe_sections[adjacent_face_indices[0]];
            let torus_section = StripSection::new(
                &pipe_section_1, &pipe_section_2, &orthogonal_pipe_section, self.bloated
            );
            torus_section.as_mesh(&config.strips)
        }).collect::<Vec<_>>()
    }

    pub fn strips_as_mesh(&self, w: f64, config: &Config) -> Mesh {
        Mesh::merge(self.strips_as_meshes(w, config))
    }

    pub fn as_mesh(&self, w: f64, config: &Config) -> ColoredMeshCollection {
        let strip_color = Color { red: 255, green: 255, blue: 255 };
        let ring_color = Color { red: 150, green: 150, blue: 150 };
        let twister_colors = vec![
            Color { red: 255, green: 0, blue: 238 },
            Color { red: 25, green: 25, blue: 255 },
        ];
        ColoredMeshCollection {
            meshes: vec![
                (self.rings_as_mesh(w, config), ring_color),
                (self.strips_as_mesh(w, config), strip_color),
                (self.twister_orbit_as_mesh(w, 0, config), twister_colors[0]),
                (self.twister_orbit_as_mesh(w, 1, config), twister_colors[1]),
            ]
        }
    }
}