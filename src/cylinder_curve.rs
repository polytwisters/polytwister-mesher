use core::f64;

use crate::{cylinder::Cylinder, pipe_section::{self, PipeSection}, polyline::Polyline};
use nalgebra::{Affine3, Point3};
use crate::utils::lerp;


enum CCurveKind {
    WrappedLoop(bool), // branch
    SideLoop(f64, f64), // theta1, theta2
    Plane(f64), // z
}


/// A CCurve is a closed curve which is one connected component of the intersection of two pipe
/// sections.
/// 
/// The curve is encoded as follows. Let C = {(x, y, z) : x^2 + y^2 = 1} be a "base cylinder." Given
/// a second cylinder B and an affine transformation A, the intersection is A*intersect(B, C). 
/// This intersection is zero, or one, or two closed curves. A CCurve is one connected component of
/// that intersection.
pub struct CCurve {
    kind: CCurveKind,
    cylinder: Cylinder,
    transform: Affine3<f64>,
}

impl CCurve {
    pub fn wrapped_loop(cylinder: Cylinder, branch: bool) -> Self {
        Self {
            kind: CCurveKind::WrappedLoop(branch),
            cylinder,
            transform: Affine3::identity(),
        } 
    }

    pub fn side_loop(cylinder: Cylinder, theta1: f64, theta2: f64) -> Self {
        Self {
            kind: CCurveKind::SideLoop(theta1, theta2),
            cylinder,
            transform: Affine3::identity(),
        } 
    }

    pub fn plane(cylinder: Cylinder, z: f64) -> Self {
        Self {
            kind: CCurveKind::Plane(z),
            cylinder,
            // NOTE: in practice "transform" is not used with CCurveKind::Plane, but it is
            // supported.
            transform: Affine3::identity(),
        } 
    }

    pub fn transform(self, transform: &Affine3<f64>) -> Self {
        Self {
            kind: self.kind,
            cylinder: self.cylinder,
            // Order doesn't really matter here since in practice self.transform is identity.
            transform: transform * self.transform
        }
    }

    /// Find the point on the curve parametrized by t, ranging form 0 to 1.
    pub fn at(&self, t: f64) -> Point3<f64> {
        let untransformed_point = match self.kind {
            CCurveKind::Plane(z) => {
                let theta = t * f64::consts::TAU;
                self.cylinder.intersect_z_plane_parametrized(z, theta)
            },
            CCurveKind::WrappedLoop(branch) => {
                let theta = t * f64::consts::TAU;
                let (p1, p2) = self.cylinder.intersect_z_line_theta(theta);
                if branch { p1 } else { p2 }
            }
            CCurveKind::SideLoop(theta1, theta2) => {
                if t < 0.5 {
                    let t2 = t * 2.0;
                    let theta = lerp(theta1, theta2, t2);
                    let (p1, _) = self.cylinder.intersect_z_line_theta(theta);
                    p1
                } else {
                    let t2 = (1.0 - t) * 2.0;
                    let theta = lerp(theta1, theta2, t2);
                    let (_, p2) = self.cylinder.intersect_z_line_theta(theta);
                    p2
                }
            },
        };
        self.transform.transform_point(&untransformed_point)
    }

    /// Given a point p in 3D space, find a value of t so that curve.at(t) is close to p.
    pub fn to_t(&self, p: &Point3<f64>) {
        let p2 = self.transform.inverse_transform_point(p);
        let untransformed_point = match self.kind {
            CCurveKind::Plane(z) => {
                
            },
            CCurveKind::WrappedLoop(branch) => {
            }
            CCurveKind::SideLoop(theta1, theta2) => {
            },
        };
    }

    pub fn discretize(&self, resolution: usize) -> Polyline {
        let points = (0..resolution).map(|i| {
            let t = i as f64 / resolution as f64;
            self.at(t)
        }).collect::<Vec<_>>();
        Polyline { points }
    }
}