use core::f64;
use std::mem::Discriminant;

use crate::{config::{CylinderMeshConfig, TorusMeshConfig}, cylinder::Cylinder, mesh::Mesh, pipe_section::{self, PipeSection}, polyline::Polyline, utils::{bisection_search, linspace, sort2}};
use nalgebra::{Affine3, Point3, Point2};
use crate::utils::{lerp, lerp_inverse};


#[derive(Clone, Copy, Debug)]
pub enum CCurveKind {
    WrappedLoop(bool), // branch
    SideLoop(f64, f64), // theta1, theta2
    Plane(f64), // z
}


/// A CCurve is a closed curve which is one connected component of the intersection of two pipe
/// sections.
/// 
/// The curve is encoded as follows. Let C = {(x, y, z) : x^2 + y^2 = 1} be a "base cylinder." Given
/// a second cylinder B and an invertible affine transformation A, the intersection is
/// A*intersect(B, C). This intersection is zero, or one, or two closed curves. A CCurve is one
/// connected component of that intersection.
/// 
/// Alternatively the curve is the intersection with a plane with a given z-coordinate.
#[derive(Clone, Copy, Debug)]
pub struct CCurve {
    pub kind: CCurveKind,
    cylinder: Cylinder,
    transform: Affine3<f64>,
}

fn warp_semicircle(x: f64) -> f64 {
    (1.0 - (x * f64::consts::PI).cos()) / 2.0
}

fn unwarp_semicircle(x: f64) -> f64 {
    (1.0 - x * 2.0).acos() * f64::consts::FRAC_1_PI
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
        let t2 = t.rem_euclid(1.0);
        let untransformed_point = match self.kind {
            CCurveKind::Plane(z) => {
                let theta = t2 * f64::consts::TAU;
                self.cylinder.intersect_z_plane_parametrized(z, theta)
            },
            CCurveKind::WrappedLoop(branch) => {
                let theta = t2 * f64::consts::TAU;
                let (p1, p2) = self.cylinder.intersect_z_line_theta(theta);
                if branch { p1 } else { p2 }
            }
            CCurveKind::SideLoop(theta1, theta2) => {
                if t2 < 0.5 {
                    // Map [0, 0.5] -> [0, 1].
                    let tmp = t2 * 2.0;
                    let tmp = warp_semicircle(tmp);
                    let theta = lerp(theta1, theta2, tmp);
                    let (p1, _) = self.cylinder.intersect_z_line_theta(theta);
                    p1
                } else {
                    // Map [0.5, 1] -> [1, 0].
                    let tmp = t2 * 2.0 - 1.0;
                    let tmp = 1.0 - tmp;
                    let tmp = warp_semicircle(tmp);
                    let theta = lerp(theta1, theta2, tmp);
                    let (_, p2) = self.cylinder.intersect_z_line_theta(theta);
                    p2
                }
            },
        };
        self.transform.transform_point(&untransformed_point)
    }

    /// Return true if the given 3D point is on the curve.
    pub fn contains(&self, p: &Point3<f64>) -> bool {
        let tolerance = 1e-5;
        let p2 = self.transform.inverse_transform_point(p);
        let x = p2.x;
        let y = p2.y;
        let z = p2.z;
        let theta = y.atan2(x).rem_euclid(f64::consts::TAU);
        let (_, z1, z2) = self.cylinder.intersect_z_line_core(&Point2::new(x, y));
        match self.kind {
            CCurveKind::Plane(plane_z) => {
                self.cylinder.scalar_field(&p2).abs() < tolerance
                && (z - plane_z).abs() < tolerance
            },
            CCurveKind::WrappedLoop(branch) => {
                (x.hypot(y) - 1.0).abs() < tolerance && if branch {
                    (z - z1).abs() < tolerance
                } else {
                    (z - z2).abs() < tolerance
                }
            }
            CCurveKind::SideLoop(theta1, theta2) => {
                let unwrapped_theta = if theta <= theta1 - tolerance {
                    theta + f64::consts::TAU
                } else {
                    theta
                };

                (x.hypot(y) - 1.0).abs() < tolerance
                && theta1 - tolerance <= unwrapped_theta
                && unwrapped_theta <= theta2 + tolerance
                && (
                    (z - z1).abs() < tolerance
                    || (z - z2).abs() < tolerance
                )
            }
        }
    }

    /// Given a point p in 3D space, find a value of t so that curve.at(t) is close to p.
    pub fn to_t(&self, p: &Point3<f64>) -> f64 {
        let p2 = self.transform.inverse_transform_point(p);
        let x = p2.x;
        let y = p2.y;
        let z = p2.z;
        let theta = y.atan2(x).rem_euclid(f64::consts::TAU);
        let (_, z1, z2) = self.cylinder.intersect_z_line_core(&Point2::new(x, y));
        match self.kind {
            CCurveKind::Plane(z) => {
                self.cylinder.z_plane_ellipse_theta(z, &p2.xy()) / f64::consts::TAU
            },
            CCurveKind::WrappedLoop(branch) => {
                theta / f64::consts::TAU
            }
            CCurveKind::SideLoop(theta1, theta2) => {
                let theta_unwrapped = if theta < theta1 {
                    theta + f64::consts::TAU
                } else {
                    theta
                };
                if (z - z1).abs() < (z - z2).abs() {
                    let tmp = lerp_inverse(theta1, theta2, theta_unwrapped);
                    let tmp = tmp.clamp(0.0, 1.0);
                    let tmp = unwarp_semicircle(tmp);
                    let t = tmp / 2.0;
                    t
                } else {
                    let tmp = lerp_inverse(theta1, theta2, theta_unwrapped);
                    let tmp = tmp.clamp(0.0, 1.0);
                    let tmp = unwarp_semicircle(tmp);
                    let tmp = 1.0 - tmp;
                    let tmp = tmp / 2.0;
                    let t = tmp + 0.5;
                    // Make sure 1.0 is wrapped back to 0.0.
                    t.rem_euclid(1.0)
                }
            }
        }
    }

    pub fn discretize_segment(&self, t1: f64, t2: f64, resolution: usize) -> Polyline {
        let points = (0..resolution).map(|i| {
            let t = lerp(t1, t2, i as f64 / (resolution as f64 - 1.0));
            self.at(t)
        }).collect::<Vec<_>>();
        Polyline { points, closed: false }
    }

    pub fn discretize_full(&self, resolution: usize) -> Polyline {
        let points = (0..resolution).map(|i| {
            let t = i as f64 / resolution as f64;
            self.at(t)
        }).collect::<Vec<_>>();
        Polyline { points, closed: true }
    }
}

#[derive(Clone, Debug)]
pub struct CylinderIntersection {
    pub ccurves: Vec<CCurve>
}

impl CylinderIntersection {
    pub fn empty() -> Self {
        Self { ccurves: vec![] }
    }

    pub fn as_mesh(&self, config: &TorusMeshConfig) -> Mesh {
        Mesh::merge(self.ccurves.iter().map(|ccurve|
            ccurve.discretize_full(config.linear_segments)
                .as_mesh(config.thickness, config.radial_segments)
        ).collect::<_>())
    }

    pub fn transform(&self, transform: &Affine3<f64>) -> Self {
        Self {
            ccurves: self.ccurves.iter().map(|ccurve|
                ccurve.transform(&transform)
            ).collect::<_>()
        } 
    }
}

#[cfg(test)]
mod test {
    use approx::*;
    use crate::{cylinder::{self, Cylinder}, cylinder_curve::{CCurveKind, unwarp_semicircle, warp_semicircle}, mesh::MeshLike, pipe_section::PipeSection, utils::linspace};
    use crate::mesh::Mesh;
    use crate::config::CylinderMeshConfig;
    use std::path::PathBuf;

    fn example_cylinder() -> Cylinder {
        Cylinder {
            m11: 1.2,
            m12: -0.5,
            m13: 0.4,
            m14: -0.1,
            m21: 0.5,
            m22: -1.2,
            m23: 0.2,
            m24: 0.05,
        }
    }

    #[test]
    fn test_to_t_plane() {
        let curve = example_cylinder().intersect_z_plane(2.0);
        let t = 0.34;
        let p = curve.at(t);
        assert!(curve.contains(&p));
        assert_abs_diff_eq!(curve.to_t(&p), t);
    }

    #[test]
    fn test_to_t_plane_2() {
        let curve = example_cylinder().intersect_z_plane(-2.0);
        let t = 0.34;
        let p = curve.at(t);
        assert!(curve.contains(&p));
        assert_abs_diff_eq!(curve.to_t(&p), t);
    }

    #[test]
    fn test_to_t_wrapped_loop() {
        let cylinder_1 = example_cylinder();
        let cylinder_2 = Cylinder::base();
        let intersection = cylinder_1.intersect_cylinder(&cylinder_2);
        let curve = intersection.ccurves[0];
        assert!(matches!(curve.kind, CCurveKind::WrappedLoop(_)));
        for t in [0.0, 0.14, 0.5, 0.99] {
            let p = curve.at(t);
            assert!(curve.contains(&p));
            assert_abs_diff_eq!(curve.to_t(&p), t);
        }
    }

    #[test]
    fn test_to_t_side_loop() {
        let cylinder_1 = Cylinder::base();
        let cylinder_2 = example_cylinder();
        let intersection = cylinder_1.intersect_cylinder(&cylinder_2);

        let curve = intersection.ccurves[0];
        assert!(matches!(curve.kind, CCurveKind::SideLoop(_, _)));
        for t in [0.0, 0.023, 0.5, 0.58] {
            let p = curve.at(t);
            assert!(curve.contains(&p));
            assert_abs_diff_eq!(curve.to_t(&p), t);
        }
    }

    #[test]
    fn test_warp_unwarp() {
        let t = 0.34;
        let w = warp_semicircle(t);
        assert!(0.0 < w && w < 1.0);
        assert_abs_diff_eq!(unwarp_semicircle(w), t);
    }

    #[test]
    fn test_unwarp_warp() {
        let t = 0.94;
        let w = unwarp_semicircle(t);
        assert!(0.0 < w && w < 1.0);
        assert_abs_diff_eq!(warp_semicircle(w), t);
    }
}