use core::f64;
use crate::{c2::C2, config::RingMeshConfig, mesh::Mesh, pipe_section::{self, PipeSection}, ring::RingSection};
use crate::config::CylinderMeshConfig;
use crate::marching_squares::{Isosurface, meshify, Grid, GridAxis};
use na::{Vector4, Point3, Vector3};

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

    pub fn twisters_as_mesh(&self, w: f64) -> Mesh {
        let config = CylinderMeshConfig {
            linear_segments: 50,
            radial_segments: 50,
            half_length: 2.0
        };
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
        Mesh::merge(meshes)
    }

    pub fn rings_as_mesh(&self, w: f64) -> Mesh {
        let config = RingMeshConfig {
            latitudes: 16,
            longitudes: 16,
            radius: 0.1,
        };
        let mut meshes = vec![];
        for ring in self.rings.iter() {
            let ring_section = RingSection {
                a: ring.vec.x.re,
                b: ring.vec.x.im,
                c: ring.vec.y.re,
                d: ring.vec.y.im,
                w
            };
            let mesh = ring_section.as_mesh(&config);
            meshes.push(mesh);
        }
        Mesh::merge(meshes)
    }


    pub fn as_mesh(&self, w: f64) -> Mesh {
        Mesh::merge(vec![
            self.rings_as_mesh(w),
            self.twisters_as_mesh(w)
        ])
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
            C2::from_parts(1.0, 0.0, 0.2, 0.3),
            C2::from_parts(-0.3, 0.4, 0.4, 0.8),
            C2::from_parts(0.7, -0.2, 0.5, 0.4),
        ]);
        let w = 0.1;
        polytwister.as_mesh(w).write_ply_file(&PathBuf::from("convex.ply"));
    }
}