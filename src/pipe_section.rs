extern crate nalgebra as na;
use na::{Point3, Vector3, Vector4};
use crate::cylinder::{Cylinder};
use crate::config::{CylinderMeshConfig, TorusMeshConfig};
use crate::mesh::{Mesh};
use crate::pipe_section;
use crate::polytwister::{FillingRegion, RegionMode};
use crate::utils::{squared};

/**
 * A 3D cross section of a pipe. (PipeCrossSection felt too long.)
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
}


/// Cross section of a twister. Currently only for convex twisters.
pub struct TwisterSection {
    pub pipe_section: PipeSection,
    pub orthogonal_pipe_section: PipeSection,
    pub neighboring_pipe_sections: Vec<PipeSection>,
    pub filling: Vec<FillingRegion>,
}

impl TwisterSection {
    pub fn contains(&self, point: &Point3<f64>) -> bool {
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

    pub fn as_mesh(&self, config: &CylinderMeshConfig) -> Mesh {
        let mut mesh = self.pipe_section.as_mesh_partial(
            |p| self.contains(p),
            &config
        );
        mesh
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
        let polylines = if self.pipe_section_1.is_plane() {
            if self.pipe_section_2.is_plane() {
                // Both are planes. Intersection is empty.
                vec![]
            } else {
                // Pipe section 1 is plane but pipe section 2 is not.
                if let Some(z) = self.pipe_section_1.plane_z() {
                    self.pipe_section_2.as_cylinder().intersect_z_planes(z, config.linear_segments)
                } else {
                    vec![]
                }
            }
        } else {
            if self.pipe_section_2.is_plane() {
                // Pipe section 2 is plane but pipe section 1 is not.
                if let Some(z) = self.pipe_section_2.plane_z() {
                    self.pipe_section_1.as_cylinder().intersect_z_planes(z, config.linear_segments)
                } else {
                    vec![]
                }
            } else {
                // Neither are planes.
                self.pipe_section_1.as_cylinder().intersect_cylinder(
                    &self.pipe_section_2.as_cylinder(), config.linear_segments
                )
            }
        };

        let meshes = polylines.iter().map(|polyline|
            polyline.as_mesh_partial(|p| {
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