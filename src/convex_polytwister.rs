use core::f64;
use serde::Deserialize;
use crate::{c2::C2, mesh::Mesh, pipe_section::{self, PipeSection}, polytwister::Polytwister, ring::RingSection, utils::linspace};
use crate::config::{CylinderMeshConfig, RingMeshConfig, TorusMeshConfig};
use crate::marching_squares::{Isosurface, meshify, Grid, GridAxis};
use crate::cylinder_curve::CCurve;
use crate::config::Config;
use na::{Vector4, Point3, Vector3, Point4};

#[derive(Deserialize, Clone)]
#[serde(rename_all="camelCase")]
pub struct ConvexPolytwisterSpec {
    logs: Vec<Vector4<f64>>,
}

pub struct ConvexPolytwister {
    logs: Vec<C2>,
    rings: Vec<C2>,
}

impl ConvexPolytwisterSpec {
    fn as_convex_polytwister(&self) -> ConvexPolytwister {
        ConvexPolytwister::new(self.logs.iter().map(|log| C2::from_vector4(&log)).collect::<_>())
    }
}

impl ConvexPolytwister {
    fn logs_contains(logs: &Vec<C2>, fiber: &C2, skip_1: usize, skip_2: usize, skip_3: usize) -> bool {
        for (index, log) in logs.iter().enumerate() {
            if index == skip_1 || index == skip_2 || index == skip_3 {
                continue;
            }
            if log.inner_abs(&fiber) > 1.0 {
                return false;
            }
        }
        true
    }

    fn deduplicate_rings(rings: &Vec<C2>) -> Vec<C2> {
        let epsilon = 1e-5;
        let mut result = vec![];
        for ring in rings.iter() {
            if result.iter().all(|ring2| ring.similarity(ring2) < 1.0 - epsilon) {
                result.push(ring.clone());
            }
        }
        result
    }

    fn compute_rings(logs: &Vec<C2>) -> Vec<C2> {
        let num_logs = logs.len();
        let mut result = vec![];
        for i in 0..num_logs {
            for j in (i + 1)..num_logs {
                for k in (j + 1)..num_logs {
                    let pipe1 = logs[i];
                    let pipe2 = logs[j];
                    let pipe3 = logs[k];
                    if let Some((ring1, ring2)) = C2::intersect_pipes(&pipe1, &pipe2, &pipe3) {
                        if Self::logs_contains(&logs, &ring1, i, j, k) {
                            result.push(ring1);
                        }
                        if Self::logs_contains(&logs, &ring2, i, j, k) {
                            result.push(ring2);
                        }
                    }
                }
            }
        }
        Self::deduplicate_rings(&result)
    }

    pub fn new(logs: Vec<C2>) -> Self {
        let rings = Self::compute_rings(&logs);
        ConvexPolytwister { logs, rings }
    }

    fn section_contains(&self, w: f64, point: &Point3<f64>) -> bool {
        for log in self.logs.iter() {
            let log_section = PipeSection::from_vector4(&log.to_vector4(), w);
            if !log_section.interior_contains(&point) {
                return false;
            }
        }
        true
    }

    fn section_contains_skip2(&self, w: f64, point: &Point3<f64>, skip_1: usize, skip_2: usize) -> bool {
        for (i, log) in self.logs.iter().enumerate() {
            if i == skip_1 || i == skip_2 {
                continue;
            }
            let log_section = PipeSection::from_vector4(&log.to_vector4(), w);
            if !log_section.interior_contains(&point) {
                return false;
            }
        }
        true
    }

    fn mesh_strip(&self, ccurve: CCurve, w: f64, config: &TorusMeshConfig, skip_1: usize, skip_2: usize) -> Mesh {
        let mut points = vec![];
        for ring in self.rings.iter() {
            let section = RingSection::from_vector4(&ring.to_vector4(), w);
            section.add_points_to_vec(&mut points);
        }
        let mut t_values = points.iter().filter_map(|point|
            if ccurve.contains(&point) {
                Some(ccurve.to_t(&point))
            } else {
                None
            }
        ).collect::<Vec<_>>();
        t_values.sort_by(f64::total_cmp);
        let mut meshes = vec![];
        for i in 0..t_values.len() {
            let t1 = t_values[i];
            let t2 = if i == t_values.len() - 1 { t_values[0] + 1.0 } else { t_values[i + 1] };
            let t_test = (t1 + t2) / 2.0;
            let p_test = ccurve.at(t_test);
            let contains = self.section_contains_skip2(w, &p_test, skip_1, skip_2);
            dbg!(i, contains);
            if contains {
                let polyline = ccurve.discretize_segment(t1, t2, config.linear_segments);
                let mesh = polyline.as_mesh(config.radius, config.radial_segments);
                meshes.push(mesh);
            }
        }
        Mesh::merge(meshes)
    }
}

impl Polytwister for ConvexPolytwister {
    fn twister_orbit_as_meshes(&self, w: f64, orbit: u8, config: &Config) -> Vec<Mesh> {
        let config = config.twisters;
        let mut meshes = vec![];
        for (i, pipe) in self.logs.iter().enumerate() {
            let pipe = pipe.rotate_real_b();
            let pipe_section = PipeSection::from_vector4(&pipe.to_vector4(), w);
            let grid = Grid {
                u_axis: GridAxis::Linear(config.linear_segments, -config.half_length, config.half_length),
                v_axis: GridAxis::Circular(config.radial_segments, f64::consts::TAU),
            };
            let surface = TwisterSection {
                pipe_section,
                log_sections: self.logs.iter().enumerate().filter_map(|(j, log)|
                    if (i == j) {
                        None
                    } else {
                        Some(PipeSection::from_vector4(&log.to_vector4(), w))
                    }
                ).collect::<_>()
            };
            let mesh = meshify(&surface, &grid);
            meshes.push(mesh);
        }
        meshes
    }

    fn strips_as_meshes(&self, w: f64, config: &Config) -> Vec<Mesh> {
        let config = config.strips;
        let mut meshes = vec![];
        for i in 0..self.logs.len() {
            for j in (i + 1)..self.logs.len() {
                let pipe1 = PipeSection::from_vector4(&self.logs[i].to_vector4(), w);
                let pipe2 = PipeSection::from_vector4(&self.logs[j].to_vector4(), w);
                let torus_section = pipe1.intersect(&pipe2);
                for ccurve in torus_section.ccurves {
                    let mesh = self.mesh_strip(ccurve, w, &config, i, j);
                    meshes.push(mesh);
                }
            }
        }
        meshes
    }

    fn rings_as_meshes(&self, w: f64, config: &Config) -> Vec<Mesh> {
        let mut meshes = vec![];
        for ring in self.rings.iter() {
            let ring_section = RingSection::from_vector4(&ring.to_vector4(), w);
            let mesh = ring_section.as_mesh(&config.rings);
            meshes.push(mesh);
        }
        meshes
    }
}

struct TwisterSection {
    pub pipe_section: PipeSection,
    pub log_sections: Vec<PipeSection>
}

impl Isosurface for TwisterSection {
    fn contains_point(&self, p: &Point3<f64>) -> bool {
        for log_section in self.log_sections.iter() {
            if !log_section.interior_contains(p) {
                return false;
            }
        }
        return true;
    }
    fn surface_coords_to_point(&self, u: f64, theta: f64) -> Point3<f64> {
        self.pipe_section.as_cylinder().surface_coords_to_cartesian(u, theta)
    }
    fn surface_coords_to_normal(&self, u: f64, theta: f64) -> Vector3<f64> {
        self.pipe_section.as_cylinder().scalar_field_gradient(
            &self.pipe_section.as_cylinder().surface_coords_to_cartesian(u, theta)
        )
    }
}

#[cfg(test)]
mod test {
    use std::path::PathBuf;
    use crate::mesh::MeshLike;
    use super::*;

    #[test]
    fn test_convex_polytwister() {
        let polytwister = ConvexPolytwister::new(vec![
            C2::from_parts(0.63, -0.42, -0.78, 0.04),
            C2::from_parts(-0.42, 0.69, 0.41, 0.32),
            C2::from_parts(-0.04, 0.54, -0.57, 0.51),
            C2::from_parts(-0.78, -0.49, -0.15, 0.10),
            C2::from_parts(0.20, 0.84, -0.51, 0.05),
            C2::from_parts(-0.27, -0.44, -0.62, -0.73)
        ]);
        let w = 0.3;
        polytwister.as_colored_mesh(w, &Config::default()).write_ply_file(&PathBuf::from("convex.ply"));
    }
}