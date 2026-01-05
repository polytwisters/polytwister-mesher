extern crate nalgebra as na;
use core::f64;

use na::{Point3, Vector3, Vector4};
use crate::cylinder::{Cylinder};
use crate::config::{CylinderMeshConfig, TorusMeshConfig};
use crate::cylinder_curve::CCurve;
use crate::mesh::{Mesh};
use crate::pipe_section;
use crate::polytwister::{FillingRegion, RegionMode};
use crate::utils::{squared};
use crate::marching_squares::{Isosurface, Grid, GridAxis, meshify};

/**
 * A 3D cross section of a pipe. May be an affine-transformed cylinder, a pair of planes, or empty.
 */
#[derive(Clone, Copy, Debug)]
pub struct PipeSection {
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub d: f64,
    pub w: f64,
}


impl PipeSection {
    pub fn new(a: f64, b: f64, c: f64, d: f64, w: f64) -> Self {
        PipeSection { a, b, c, d, w }
    }

    pub fn from_vector4(vector: &Vector4<f64>, w: f64) -> Self {
        PipeSection {
            a: vector.x,
            b: vector.y,
            c: vector.z,
            d: vector.w,
            w
        }
    }

    /**
     * Evaluate the scalar field associated with this pipe cross section:
     * 
     *     F(x, y, z) = (ax + by + cz + dw)^2 + (bx - ay + dz - cw)^2 - 1
     * 
     * The pipe cross section is given by the isosurface F(x, y, z) = 0. For other points, if
     * F(x, y, z) < 0 then the point is inside the pipe, F(x, y, z) > 0 is outside the pipe, and
     * F(x, y, z) = -1 is on the pipe's symmetry axis (plane of symmetry if a = b = 0).
     */
    pub fn scalar_field(&self, point: &Point3<f64>) -> f64 {
        squared(self.a * point.x + self.b * point.y + self.c * point.z + self.d * self.w)
        + squared(self.b * point.x - self.a * point.y + self.d * point.z - self.c * self.w)
        - 1.0
    }

    pub fn contains(&self, point: &Point3<f64>) -> bool {
        self.scalar_field(point) < 0.0
    }

    pub fn as_cylinder(&self) -> Cylinder {
        Cylinder {
            m11: self.a,
            m12: self.b,
            m13: self.c,
            m14: self.d * self.w,
            m21: self.b,
            m22: -self.a,
            m23: self.d,
            m24: -self.c * self.w,
        }
    }

    pub fn is_plane(&self) -> bool {
        self.a == 0.0 && self.b == 0.0
    }

    pub fn plane_z(&self) -> Option<f64> {
        let tmp = 1.0 / squared(self.c) - squared(self.w);
        if tmp <= 0.0 {
            None
        } else {
            Some(tmp.sqrt())
        }
    }

    pub fn as_mesh(&self, config: &CylinderMeshConfig) -> Mesh {
        if self.d != 0.0 {
            panic!("PipeSection::as_mesh does not yet work with d != 0");
        }

        if self.is_plane() {
            if let Some(z) = self.plane_z() {
                Mesh::merge(vec![
                    Mesh::plane(z, config.half_length, config.linear_segments),
                    Mesh::plane(-z, config.half_length, config.linear_segments),
                ])
            } else {
                Mesh::empty()
            }
        } else {
            self.as_cylinder().as_mesh(config)
        }
    }

    pub fn as_mesh_partial<F: Fn(&Point3<f64>) -> bool>(
        &self,
        predicate: F,
        config: &CylinderMeshConfig
    ) -> Mesh {
        if self.d != 0.0 {
            panic!("PipeSection::as_mesh does not yet work with d != 0");
        }

        if self.is_plane() {
            if let Some(z) = self.plane_z() {
                Mesh::merge(vec![
                    Mesh::partial_plane(&predicate, z, config.half_length, config.linear_segments),
                    Mesh::partial_plane(&predicate, -z, config.half_length, config.linear_segments),
                ])
            } else {
                Mesh::empty()
            }
        } else {
            self.as_cylinder().as_mesh_partial(&predicate, config)
        }
    }

    fn intersect(&self, other: &PipeSection) -> Vec<CCurve> {
        if self.is_plane() {
            if other.is_plane() {
                // Both are planes. Intersection is empty.
                vec![]
            } else {
                // Pipe section 1 is plane, pipe section 2 is cylinder.
                if let Some(z) = self.plane_z() {
                    other.as_cylinder().intersect_z_planes(z)
                } else {
                    // Pipe section 1 is empty.
                    vec![]
                }
            }
        } else {
            if other.is_plane() {
                if let Some(z) = other.plane_z() {
                    // Pipe section 1 is cylinder, pipe section 2 is plane.
                    self.as_cylinder().intersect_z_planes(z)
                } else {
                    // Pipe section 2 is empty.
                    vec![]
                }
            } else {
                // Both are cylinders.
                self.as_cylinder().intersect_cylinder(&other.as_cylinder())
            }
        }
    }
}


/// Cross section of a twister.
pub struct TwisterSection {
    pub pipe_section: PipeSection,
    pub filling_info: TwisterFillingInfo,
}

/// All the geometric information needed to fill a twister, without the twister pipe section itself.
#[derive(Clone)]
pub struct TwisterFillingInfo {
    pub orthogonal_pipe_section: PipeSection,
    pub neighboring_pipe_sections: Vec<PipeSection>,
    pub filling: Vec<FillingRegion>,
}

pub struct TwisterCylindricalIsosurface {
    pub pipe_section: PipeSection,
    pub filling_info: TwisterFillingInfo,
}

pub struct TwisterPlanarIsosurface {
    pub z: f64,
    pub filling_info: TwisterFillingInfo,
}

impl TwisterFillingInfo {
    pub fn contains_point(&self, point: &Point3<f64>) -> bool {
        let inner = self.orthogonal_pipe_section.contains(point);
        let outer = !inner;
        let order: u32 = self.neighboring_pipe_sections.iter().map(|ps|
            if ps.contains(point) { 1 } else { 0 }
        ).sum();
        self.filling.iter().any(|region| {
            region.order == order
            && match region.mode {
                RegionMode::Inner => inner,
                RegionMode::Outer => outer,
                RegionMode::Both => true
            }
        })
    }
}

impl Isosurface for TwisterCylindricalIsosurface {
    fn contains_point(&self, p: &Point3<f64>) -> bool {
        self.filling_info.contains_point(p)
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

impl Isosurface for TwisterPlanarIsosurface {
    fn contains_point(&self, p: &Point3<f64>) -> bool {
        self.filling_info.contains_point(p)
    }
    fn surface_coords_to_point(&self, u: f64, v: f64) -> Point3<f64> {
        Point3::new(u, v, self.z)
    }
    fn surface_coords_to_normal(&self, u: f64, v: f64) -> Vector3<f64> {
        Vector3::z()
    }
}

impl TwisterSection {
    pub fn new(
        pipe_section: PipeSection,
        orthogonal_pipe_section: PipeSection,
        neighboring_pipe_sections: Vec<PipeSection>,
        filling: Vec<FillingRegion>,
    ) -> Self {
        Self {
            pipe_section,
            filling_info: TwisterFillingInfo {
                orthogonal_pipe_section,
                neighboring_pipe_sections,
                filling
            }
        }
    }

    pub fn as_mesh(&self, config: &CylinderMeshConfig) -> Mesh {
        let extent = 2.0;
        let segments = 30;
        if self.pipe_section.is_plane() {
            if let Some(z) = self.pipe_section.plane_z() {
                let grid = Grid {
                    u_axis: GridAxis::Linear(segments, -extent, extent),
                    v_axis: GridAxis::Linear(segments, -extent, extent),
                };
                let plane_1 = TwisterPlanarIsosurface {
                    z: z,
                    filling_info: self.filling_info.clone(),
                };
                let plane_2 = TwisterPlanarIsosurface {
                    z: -z,
                    filling_info: self.filling_info.clone(),
                };
                Mesh::merge(vec![
                    meshify(&plane_1, &grid),
                    meshify(&plane_2, &grid),
                ])
            } else {
                Mesh::empty()
            }
        } else {
            let grid = Grid {
                u_axis: GridAxis::Linear(segments, -extent, extent),
                v_axis: GridAxis::Circular(segments, f64::consts::TAU),
            };
            let surface = TwisterCylindricalIsosurface {
                pipe_section: self.pipe_section,
                filling_info: self.filling_info.clone(),
            };
            meshify(&surface, &grid)
        }
    }
}


#[derive(Clone, Copy, Debug)]
pub struct StripSection {
    pub pipe_section_1: PipeSection,
    pub pipe_section_2: PipeSection,
    pub orthogonal_pipe_section: PipeSection,
    pub bloated: bool,
}

impl StripSection {
    pub fn new(
        pipe_section_1: &PipeSection,
        pipe_section_2: &PipeSection,
        orthogonal_pipe_section: &PipeSection,
        bloated: bool,
    ) -> Self {
        Self {
            pipe_section_1: pipe_section_1.clone(),
            pipe_section_2: pipe_section_2.clone(),
            orthogonal_pipe_section: orthogonal_pipe_section.clone(),
            bloated,
        }
    }

    pub fn as_mesh(&self, config: &TorusMeshConfig) -> Mesh {
        let ccurves = self.pipe_section_1.intersect(&self.pipe_section_2);

        let meshes = ccurves.iter().map(|ccurve|
            ccurve
                .discretize(config.linear_segments)
                .as_mesh_partial(|p| {
                    if self.bloated {
                        !self.orthogonal_pipe_section.contains(&p)
                    } else {
                        self.orthogonal_pipe_section.contains(&p)
                    }
                }, config.thickness, config.radial_segments)
        ).collect::<_>();
        Mesh::merge(meshes)
    }
}