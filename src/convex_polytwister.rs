use core::f64;
use serde::Deserialize;
use crate::{c2::C2, elements::Pipe, mesh::Mesh, pipe_section::{self, PipeSection}, polytwister::Polytwister, ring::RingSection, utils::linspace};
use crate::config::{CylinderMeshConfig, RingMeshConfig, TorusMeshConfig};
use crate::marching_squares::{Isosurface, meshify, Grid, GridAxis};
use crate::cylinder_curve::CCurve;
use crate::config::Config;
use crate::elements::{Fiber, Log};
use crate::utils::UnorderedTriples;
use na::{Vector4, Point3, Vector3, Point4};

/// A specification of a convex polytwister intended for deserialization only. It is represented as
/// a list of C^2 vectors x_i so that the polytwister is the intersection of all L(x_i), where
/// L(y) is the set of all x in C^2 such that |<x, y>| <= 1. The C^2 vectors x_i are specified as
/// Vector4s with the real and imaginary parts split.
#[derive(Deserialize, Clone)]
#[serde(rename_all="camelCase")]
pub struct ConvexPolytwisterSpec {
    logs: Vec<Vector4<f64>>,
}

pub struct ConvexPolytwister {
    logs: Vec<Log>,
    rings: Vec<Fiber>,
}

impl ConvexPolytwisterSpec {
    /// Convert this to a convex polytwister.
    pub fn to_convex_polytwister(&self) -> ConvexPolytwister {
        ConvexPolytwister::new(self.logs.iter().map(|log|
            Log::from_vector4(&log)
        ).collect::<_>())
    }
}

impl ConvexPolytwister {
    pub fn new(logs: Vec<Log>) -> Self {
        let rings = Self::compute_rings(&logs);
        ConvexPolytwister { logs, rings }
    }

    /// Given a set of logs L(p_i) and a single vector x_i, check to see whether x_i is in all logs
    /// *except* the three logs given at the indices skip_1, skip_2, or skip_3.
    fn logs_contains(logs: &Vec<Log>, fiber: &Fiber, skip_1: usize, skip_2: usize, skip_3: usize) -> bool {
        for (index, log) in logs.iter().enumerate() {
            if index == skip_1 || index == skip_2 || index == skip_3 {
                continue;
            }
            if log.contains(fiber) {
                return false;
            }
        }
        true
    }

    /// Given a set of logs L(x_i) specified using the C^2 vectors x_i, find all rings. Return each
    /// ring as a single C^2 vector.
    fn compute_rings(logs: &Vec<Log>) -> Vec<Fiber> {
        let num_logs = logs.len();
        let mut result = vec![];
        // Check all triples (i, j, k) where 0 <= i < j < k < num_logs.
        for (i, j, k) in UnorderedTriples::new(num_logs) {
            let pipe1 = logs[i].bounding_pipe();
            let pipe2 = logs[j].bounding_pipe();
            let pipe3 = logs[k].bounding_pipe();
            if let Some((ring1, ring2)) = Pipe::intersect(&pipe1, &pipe2, &pipe3) {
                if Self::logs_contains(&logs, &ring1, i, j, k) {
                    result.push(ring1);
                }
                if Self::logs_contains(&logs, &ring2, i, j, k) {
                    result.push(ring2);
                }
            }
        }
        let epsilon = 1e-5;
        Fiber::deduplicate(&result, epsilon)
    }

    fn contains(&self, fiber: &Fiber) -> bool {
        self.logs.iter().all(|log| log.contains(&fiber)) 
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
        // Sorry about code dupe with StripSection in uniform polytwisters. I was in a hurry to get this done.
        t_values.sort_by(f64::total_cmp);
        let mut meshes = vec![];
        for i in 0..t_values.len() {
            let t1 = t_values[i];
            let t2 = if i == t_values.len() - 1 { t_values[0] + 1.0 } else { t_values[i + 1] };
            let t_test = (t1 + t2) / 2.0;
            let p_test = ccurve.at(t_test);
            let contains = self.section_contains_skip2(w, &p_test, skip_1, skip_2);
            if contains {
                let polyline = ccurve.discretize_segment(t1, t2, config.linear_segments);
                let mesh = polyline.as_mesh(config.radius, config.radial_segments);
                meshes.push(mesh);
            }
        }
        Mesh::merge(meshes)
    }

    fn scale(&self, ratio: f64) -> Self {
        Self {
            logs: self.logs.iter().map(|log| log.scale(ratio)).collect::<_>(),
            rings: self.rings.iter().map(|ring| ring.scale(ratio)).collect::<_>(),
        }
    }

    fn radius(&self) -> f64 {
        // Not actually correct, but gets the job done for now
        self.rings.iter().map(|ring| ring.radius()).max_by(f64::total_cmp).unwrap_or(1.0)
    }

    pub fn normalize(&self) -> Self {
        self.scale(1.0 / self.radius())
    }
}

impl Polytwister for ConvexPolytwister {
    fn twister_orbit_as_meshes(&self, w: f64, orbit: u8, config: &Config) -> Vec<Mesh> {
        if orbit != 0 {
            return vec![];
        }

        let config = config.twisters;
        let mut meshes = vec![];
        for (i, pipe) in self.logs.iter().enumerate() {
            let pipe = pipe.vec.rotate_real_b();
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
    use super::*;

    fn dyster() -> ConvexPolytwister {
        ConvexPolytwister::new(vec![
            Log::new(C2::from_components(1.0, 0.0, 0.0, 0.0)),
            Log::new(C2::from_components(0.0, 0.0, 1.0, 0.0)),
            Log::new(C2::from_components(0.0, 1.0, 1.0, 0.0)),
        ])
    }

    #[test]
    fn test_dyster() {
        let dyster = dyster();
        let rings = dyster.rings.clone();
        assert_eq!(rings.len(), 2);
        assert!(dyster.contains(&Fiber::zero()));
        assert!(dyster.contains(&rings[0]));
        assert!(dyster.contains(&rings[1]));
    }

    #[test]
    fn test_basic() {
        let dyster = dyster();
        let meshes = dyster.as_meshes(0.3, &Config::default());
        assert!(meshes.ring_meshes.len() > 0);
        assert!(meshes.strip_meshes.len() > 0);
        assert!(meshes.twister_meshes_orbit_1.len() > 0);
        assert!(meshes.twister_meshes_orbit_2.len() == 0);
    }
}