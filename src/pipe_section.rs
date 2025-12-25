extern crate nalgebra as na;
use na::{Point3};
use crate::cylinder::{Cylinder, CylinderMeshOptions};
use crate::mesh::{Mesh};
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

    pub fn as_mesh(&self, options: &CylinderMeshOptions) -> Mesh {
        if self.d != 0.0 {
            panic!("PipeSection::as_mesh does not yet work with d != 0");
        }

        if self.is_plane() {
            if let Some(z) = self.plane_z() {
                Mesh::merge(vec![
                    Mesh::plane(z, options.half_length, options.linear_segments),
                    Mesh::plane(-z, options.half_length, options.linear_segments),
                ])
            } else {
                Mesh::empty()
            }
        } else {
            self.as_cylinder().as_mesh(options)
        }
    }

    pub fn as_mesh_partial<F: Fn(&Point3<f64>) -> bool>(
        &self,
        predicate: F,
        options: &CylinderMeshOptions
    ) -> Mesh {
        if self.d != 0.0 {
            panic!("PipeSection::as_mesh does not yet work with d != 0");
        }

        if self.is_plane() {
            if let Some(z) = self.plane_z() {
                Mesh::merge(vec![
                    Mesh::partial_plane(&predicate, z, options.half_length, options.linear_segments),
                    Mesh::partial_plane(&predicate, -z, options.half_length, options.linear_segments),
                ])
            } else {
                Mesh::empty()
            }
        } else {
            self.as_cylinder().as_mesh_partial(&predicate, options)
        }
    }
}


/// Cross section of a twister. Currently only for convex twisters.
pub struct TwisterSection {
    pub pipe_section: PipeSection,
    pub neighboring_pipe_sections: Vec<PipeSection>,
}

impl TwisterSection {
    pub fn as_mesh(&self, options: &CylinderMeshOptions) -> Mesh {
        let mut mesh = self.pipe_section.as_mesh_partial(
            |p| {
                for pipe_section_2 in self.neighboring_pipe_sections.iter() {
                    if !pipe_section_2.contains(p) {
                        return false;
                    }
                }
                true
            },
            &options
        );
        mesh
    }
}


#[derive(Clone, Copy, Debug)]
pub struct TorusSection {
    pub pipe_section_1: PipeSection,
    pub pipe_section_2: PipeSection,
}

#[derive(Clone, Copy, Debug)]
pub struct TorusMeshOptions {
    pub thickness: f64,
    pub linear_segments: usize,
    pub radial_segments: usize,
}

impl TorusSection {
    pub fn new(pipe_section_1: &PipeSection, pipe_section_2: &PipeSection) -> Self {
        Self {
            pipe_section_1: pipe_section_1.clone(),
            pipe_section_2: pipe_section_2.clone()
        }
    }

    pub fn as_mesh(&self, options: &TorusMeshOptions) -> Mesh {
        let polylines = if self.pipe_section_1.is_plane() {
            if self.pipe_section_2.is_plane() {
                // Both are planes. Intersection is empty.
                vec![]
            } else {
                // Pipe section 1 is plane but pipe section 2 is not.
                if let Some(z) = self.pipe_section_1.plane_z() {
                    self.pipe_section_2.as_cylinder().intersect_z_planes(z, options.linear_segments)
                } else {
                    vec![]
                }
            }
        } else {
            if self.pipe_section_2.is_plane() {
                // Pipe section 2 is plane but pipe section 1 is not.
                if let Some(z) = self.pipe_section_2.plane_z() {
                    self.pipe_section_1.as_cylinder().intersect_z_planes(z, options.linear_segments)
                } else {
                    vec![]
                }
            } else {
                // Neither are planes.
                self.pipe_section_1.as_cylinder().intersect_cylinder(
                    &self.pipe_section_2.as_cylinder(), options.linear_segments
                )
            }
        };
        
        let meshes = polylines.iter().map(|polyline|
            polyline.as_mesh(options.thickness, options.radial_segments)
        ).collect::<_>();
        Mesh::merge(meshes)
    }
}