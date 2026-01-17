use core::f64;
use std::mem::Discriminant;

use crate::{cylinder::Cylinder, pipe_section::{self, PipeSection}, polyline::Polyline, utils::sort2};
use nalgebra::{Affine3, Point3, Point2};
use crate::utils::{lerp, lerp_inverse};


#[derive(Clone, Copy, Debug)]
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
/// 
/// Alternatively the curve is the intersection with a plane with a given z-coordinate.
#[derive(Clone, Copy, Debug)]
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
                    // Map [0, 1] -> [0, 0.5].
                    let t3 = t2 * 2.0;
                    let theta = lerp(theta1, theta2, t3);
                    let (p1, _) = self.cylinder.intersect_z_line_theta(theta);
                    p1
                } else {
                    // Map [0.5, 1] -> [1, 0].
                    let t3 = (1.0 - t2) * 2.0;
                    let theta = lerp(theta1, theta2, t3);
                    let (_, p2) = self.cylinder.intersect_z_line_theta(theta);
                    p2
                }
            },
        };
        self.transform.transform_point(&untransformed_point)
    }

    /// Return true if the given 3D point is on the curve.
    pub fn contains(&self, p: &Point3<f64>) -> bool {
        let tolerance = 1e-10;
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
                && theta1 - tolerance <= unwrapped_theta && unwrapped_theta <= theta2 + tolerance
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
                    // Branch 1: 0 <= t < 0.5
                    lerp_inverse(theta1, theta2, theta_unwrapped).clamp(0.0, 1.0) / 2.0
                } else {
                    // Branch 2: 0.5 <= t < 1
                    let tmp = 0.5 + lerp_inverse(
                        theta2, theta1, theta_unwrapped
                    ).clamp(0.0, 1.0) / 2.0;
                    // Make sure 1.0 is wrapped back to 0.0.
                    tmp.rem_euclid(1.0)
                }
            }
        }
    }

    pub fn strip_t_interval(
        &self,
        p1: &Point3<f64>,
        p2: &Point3<f64>,
        orthogonal_pipe_section: &PipeSection,
        bloated: bool
    ) -> (f64, f64) {
        let (mut t1, mut t2) = sort2((self.to_t(&p1), self.to_t(&p2)));
        let t3 = (t1 + t2) / 2.0;
        let p3 = self.at(t3);
        let strip_contains_p3 = !orthogonal_pipe_section.interior_contains(&p3) == bloated;
        if !strip_contains_p3 {
            // Switch to complementary interval.
            (t1, t2) = (t2, t1 + 1.0);
        }
        (t1, t2)
    }

    pub fn discretize(&self, t1: f64, t2: f64, resolution: usize) -> Polyline {
        let points = (0..resolution).map(|i| {
            let t = lerp(t1, t2, i as f64 / resolution as f64);
            self.at(t)
        }).collect::<Vec<_>>();
        Polyline { points }
    }
}

#[cfg(test)]
mod test {
    use approx::*;
    use crate::{cylinder::{self, Cylinder}, cylinder_curve::CCurveKind, pipe_section::PipeSection, utils::linspace};

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
        let curves = cylinder_1.intersect_cylinder(&cylinder_2);
        let curve = curves[0];
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
        let curves = cylinder_1.intersect_cylinder(&cylinder_2);
        let curve = curves[0];
        assert!(matches!(curve.kind, CCurveKind::SideLoop(_, _)));
        for t in [0.0, 0.023, 0.5, 0.58] {
            let p = curve.at(t);
            assert!(curve.contains(&p));
            assert_abs_diff_eq!(curve.to_t(&p), t);
        }
    }

    #[test]
    fn test_intersect_2() {
        let w = 0.1;
        let pipe_section_1 = PipeSection::new(
            0.5000000000000001,
            -0.8660254037844386,
            0.5176380902050415,
            0.0,
            w,
        );
        let pipe_section_2 = PipeSection::new(
            0.5000000000000001,
            -0.8660254037844386,
            0.5176380902050415,
            0.0,
            w,
        );
        let cylinder_1 = pipe_section_1.as_cylinder();
        let cylinder_2 = pipe_section_2.as_cylinder();
        let curves = cylinder_1.intersect_cylinder(&cylinder_2);
        for curve in curves {
            for t in linspace(0.0, 1.0, 500) {
                let p = curve.at(t);
                let tolerance = 1e-5;
                assert!(pipe_section_1.boundary_contains(&p, tolerance));
                assert!(pipe_section_2.boundary_contains(&p, tolerance));
            }
        }
    }
}