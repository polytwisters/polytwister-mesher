extern crate nalgebra as na;
use core::f64;

use na::{Complex, Point3, Vector3, Vector4};
use crate::cylinder::{Cylinder};
use crate::config::{CylinderMeshConfig, TorusMeshConfig};
use crate::cylinder_curve::{CCurve, CylinderIntersection};
use crate::mesh::{Mesh};
use crate::polyline::Polyline;
use crate::{pipe_section, ring};
use crate::polytwister::{FillingRegion, RegionMode};
use crate::ring::{RingSection, RingSectionResult};
use crate::utils::{bisection_search, sort4, squared};
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

    /// Return true if the point is on the boundary or interior of the log.
    pub fn interior_contains(&self, point: &Point3<f64>) -> bool {
        self.scalar_field(point) <= 0.0
    }

    /// Return true if the point is on the boundary with the given tolerance.
    pub fn boundary_contains(&self, point: &Point3<f64>, tolerance: f64) -> bool {
        self.scalar_field(point).abs() <= tolerance
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

    pub fn intersect(&self, other: &PipeSection) -> CylinderIntersection {
        if self.is_plane() {
            if other.is_plane() {
                // Both are planes. Intersection is empty.
                CylinderIntersection::empty()
            } else {
                // Pipe section 1 is plane, pipe section 2 is cylinder.
                if let Some(z) = self.plane_z() {
                    other.as_cylinder().intersect_z_planes(z)
                } else {
                    // Pipe section 1 is empty.
                    CylinderIntersection::empty()
                }
            }
        } else {
            if other.is_plane() {
                if let Some(z) = other.plane_z() {
                    // Pipe section 1 is cylinder, pipe section 2 is plane.
                    self.as_cylinder().intersect_z_planes(z)
                } else {
                    // Pipe section 2 is empty.
                    CylinderIntersection::empty()
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
        let inner = self.orthogonal_pipe_section.interior_contains(point);
        let outer = !inner;
        let order: u32 = self.neighboring_pipe_sections.iter().map(|ps|
            if ps.interior_contains(point) { 1 } else { 0 }
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
        let extent = config.half_length;
        let linear_segments = config.linear_segments;
        let radial_segments = config.radial_segments;
        if self.pipe_section.is_plane() {
            if let Some(z) = self.pipe_section.plane_z() {
                let grid = Grid {
                    u_axis: GridAxis::Linear(linear_segments, -extent, extent),
                    v_axis: GridAxis::Linear(linear_segments, -extent, extent),
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
                u_axis: GridAxis::Linear(linear_segments, -extent, extent),
                v_axis: GridAxis::Circular(radial_segments, f64::consts::TAU),
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
    pub ring_section_1: RingSection,
    pub ring_section_2: RingSection,
    pub bloated: bool,
}

enum StripIntervals {
    None,
    All,
    OneInterval((f64, f64)),
    TwoIntervals((f64, f64), (f64, f64)),
}

impl StripSection {
    pub fn new(
        pipe_section_1: &PipeSection,
        pipe_section_2: &PipeSection,
        orthogonal_pipe_section: &PipeSection,
        ring_section_1: &RingSection,
        ring_section_2: &RingSection,
        bloated: bool,
    ) -> Self {
        Self {
            pipe_section_1: pipe_section_1.clone(),
            pipe_section_2: pipe_section_2.clone(),
            orthogonal_pipe_section: orthogonal_pipe_section.clone(),
            ring_section_1: ring_section_1.clone(),
            ring_section_2: ring_section_2.clone(),
            bloated,
        }
    }

    /// Given a point that is on the torus cross section, return whether it is also on the strip
    /// cross section by checking it against the orthogonal pipe section.
    pub fn contains_point_on_torus(&self, point: &Point3<f64>) -> bool {
        !self.orthogonal_pipe_section.interior_contains(point) == self.bloated
    }

    pub fn as_mesh(&self, config: &TorusMeshConfig) -> Mesh {
        // First get all the CCurves, curves equal to the intersection of the two pipes, and
        // therefore the cross section of the torus containing this strip.
        let intersection = self.pipe_section_1.intersect(&self.pipe_section_2);
        let ccurves = intersection.ccurves;

        // The strip cross section's endpoints are always the points which are cross sections of its
        // bounding rings. However, we do not know which ones yet. First, let's gather all the
        // candidates for endpoints of the strip cross section.
        let mut endpoints = vec![];
        let points_1 = self.ring_section_1.add_points_to_vec(&mut endpoints);
        let points_2 = self.ring_section_2.add_points_to_vec(&mut endpoints);

        let has_xy_ring = (
            matches!(self.ring_section_1.as_points(), RingSectionResult::XYCircle { radius })
            || matches!(self.ring_section_2.as_points(), RingSectionResult::XYCircle { radius })
        );

        let meshes = ccurves.iter().map(|ccurve| {
            // Some of the ring cross sections will be on the CCurve and others will not. Filter out
            // the endpoint candidates that aren't on the CCurve, and convert them to t-values in
            // the CCurve's explicit parametrization.
            let mut t_values = endpoints.iter().filter_map(|point| {
                if ccurve.contains(&point) {
                    Some(ccurve.to_t(&point))
                } else {
                    None
                }
            }).collect::<Vec<_>>();

            if has_xy_ring {
                // Stupid. Just discretize the strip with a combined grid/bisection search.
                let resolution = 1000;
                for i in 0..resolution {
                    let t1 = i as f64 / resolution as f64;
                    let t2 = (i + 1) as f64 / resolution as f64;
                    let p1 = ccurve.at(t1);
                    let p2 = ccurve.at(t2);
                    if self.contains_point_on_torus(&p1) != self.contains_point_on_torus(&p2) {
                        let t_frac = bisection_search(|t_frac| {
                            self.contains_point_on_torus(
                                &ccurve.at((i as f64 + t_frac) / resolution as f64)
                            )
                        });
                        t_values.push((i as f64 + t_frac) / resolution as f64);
                    }
                }
            }

            // t-values are in the range [0, 1], treated circularly. To reduce casework we have them
            // in increasing order.
            t_values.sort_by(f64::total_cmp);

            if t_values.len() == 0 {
                // If there are no ring cross sections on this CCurve, then either the entire CCurve
                // is part of the strip cross section, or none of it. To test this we just try one
                // point.
                let t_test = 0.25; // doesn't matter
                let p_test = ccurve.at(t_test);
                let contains = self.contains_point_on_torus(&p_test);
                if contains {
                    let polyline = ccurve.discretize_full(config.linear_segments);
                    polyline.as_mesh(config.radius, config.radial_segments)
                } else {
                    Mesh::empty()
                }
            } else {
                // Otherwise, use the t_values to divide up the CCurve into segments. Arbitrarily
                // test the midpoint of each section.
                let mut meshes = vec![];
                for i in 0..t_values.len() {
                    let t1 = t_values[i];
                    let t2 = if i == t_values.len() - 1 { t_values[0] + 1.0 } else { t_values[i + 1] };
                    let t_test = (t1 + t2) / 2.0;
                    let p_test = ccurve.at(t_test);
                    let contains = self.contains_point_on_torus(&p_test);
                    if contains {
                        let polyline = ccurve.discretize_segment(t1, t2, config.linear_segments);
                        let mesh = polyline.as_mesh(config.radius, config.radial_segments);
                        meshes.push(mesh);
                    }
                }
                Mesh::merge(meshes)
            }
        }).collect::<_>();
        Mesh::merge(meshes)
    }
}