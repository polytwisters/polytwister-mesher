use std::path::PathBuf;
use crate::config::Config;
use crate::mesh::{Mesh, MeshLike, ColoredMesh, Color};

pub trait Polytwister {
    fn rings_as_meshes(&self, w: f64, config: &Config) -> Vec<Mesh>;
    fn strips_as_meshes(&self, w: f64, config: &Config) -> Vec<Mesh>;
    fn twister_orbit_as_meshes(&self, w: f64, orbit: u8, config: &Config) -> Vec<Mesh>;

    fn as_meshes(&self, w: f64, config: &Config) -> PolytwisterMeshes {
        PolytwisterMeshes {
            ring_meshes: self.rings_as_meshes(w, config),
            strip_meshes: self.strips_as_meshes(w, config),
            twister_meshes_orbit_1: self.twister_orbit_as_meshes(w, 0, config),
            twister_meshes_orbit_2: self.twister_orbit_as_meshes(w, 1, config),
        }
    }
    fn as_colored_mesh(&self, w: f64, config: &Config) -> ColoredMesh {
        self.as_meshes(w, config).as_colored_mesh()
    }
}


pub struct PolytwisterMeshes {
    pub ring_meshes: Vec<Mesh>,
    pub strip_meshes: Vec<Mesh>,
    pub twister_meshes_orbit_1: Vec<Mesh>,
    pub twister_meshes_orbit_2: Vec<Mesh>,
}

impl PolytwisterMeshes {
    pub fn write_plys(&self, dir: &PathBuf) -> std::io::Result<()> {
        let ring_mesh = Mesh::merge(self.ring_meshes.clone());
        let strip_mesh = Mesh::merge(self.strip_meshes.clone());
        let twister_mesh_orbit_1 = Mesh::merge(self.twister_meshes_orbit_1.clone());
        let twister_mesh_orbit_2 = Mesh::merge(self.twister_meshes_orbit_2.clone());

        ring_mesh.write_ply_file_and_log(&dir.join(format!("rings.ply")), "Ring mesh")?;
        strip_mesh.write_ply_file_and_log(&dir.join(format!("strips.ply")), "Strip mesh")?;
        twister_mesh_orbit_1.write_ply_file_and_log(&dir.join(format!("twisters_1.ply")), "Twister orbit 1 mesh")?;
        twister_mesh_orbit_2.write_ply_file_and_log(&dir.join(format!("twisters_2.ply")), "Twister orbit 2 mesh")?;

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