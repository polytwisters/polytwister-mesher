use core::f64;
use std::io::Empty;

use crate::{cylinder::Cylinder, cylinder_curve::CylinderIntersection, pipe_section::{self, PipeSection}, polyline::Polyline, utils::{angle_vector, linspace}};
use crate::cylinder_curve::{CCurve};
use na::{Matrix2, Point2, Point3, Vector2, Vector3, Affine3, Matrix4};
use crate::utils::{squared, sort2, sort4, angle, unzip_circle};

/**
 * A 2D line given by {p + dt | t in R} where p and d are in R^2 and d is a unit vector.
 */
pub struct Line2D {
    pub p: Point2<f64>,
    pub d: Vector2<f64>,
}

impl Line2D {
    /**
     * Return the intersection of this line with the circle cos(theta)^2 + sin(theta)^2 = 0.
     * Returns either a tuple of two thetas (the same point twice if the line is tangent), or None
     * if the line does not intersect the circle.
     */
    pub fn intersect_unit_circle(&self) -> Option<(Point2<f64>, Point2<f64>)> {
        let px = self.p.x;
        let py = self.p.y;
        let dx = self.d.x;
        let dy = self.d.y;
        let b = 2.0 * (px * dx + py * dy);
        let c = px * px + py * py - 1.0;
        let discriminant = b * b - 4.0 * c;
        if discriminant < 0.0 {
            return None;
        };
        let tmp = discriminant.sqrt();
        let t1 = (-b + tmp) / 2.0;
        let t2 = (-b - tmp) / 2.0;
        let p1 = self.p + t1 * self.d;
        let p2 = self.p + t2 * self.d;
        Some((p1, p2))
    }
}

/**
 * A *closed* interval over angles from 0 to 2pi. We always have 0 <= x < start, but end may be
 * 2pi or greater.
 */
#[derive(Debug)]
pub struct CircularInterval {
    pub start: f64,
    pub end: f64
}

impl CircularInterval {
    fn new(start: f64, end: f64) -> Self {
        CircularInterval { start, end }
    }

    fn contains(&self, x: f64) -> bool {
        if self.end >= f64::consts::TAU {
            self.start <= x || x <= self.end - f64::consts::TAU
        } else {
            self.start <= x && x <= self.end
        }
    }
}

#[derive(Debug)]
pub enum CylinderIntersectionSolutions {
    Empty,
    All,
    OneInterval(CircularInterval),
    TwoIntervals(CircularInterval, CircularInterval),
}

impl CylinderIntersectionSolutions {
    pub fn contains(&self, theta: f64) -> bool {
        match self {
            CylinderIntersectionSolutions::Empty => false,
            CylinderIntersectionSolutions::All => true,
            CylinderIntersectionSolutions::OneInterval(interval) => interval.contains(theta),
            CylinderIntersectionSolutions::TwoIntervals(interval1, interval2) => {
                interval1.contains(theta) || interval2.contains(theta)
            },
        }
    }

    pub fn values(&self) -> Vec<f64> {
        match self {
            CylinderIntersectionSolutions::Empty => vec![],
            CylinderIntersectionSolutions::All => vec![],
            CylinderIntersectionSolutions::OneInterval(interval) => vec![interval.start, interval.end],
            CylinderIntersectionSolutions::TwoIntervals(interval1, interval2) => vec![
                interval1.start, interval1.end, interval2.start, interval2.end
            ],
        }
    }
}

impl Cylinder {

    /**
     * Given a 2D point (x, y), intersect the pipe section with the line parallel to the z-axis
     * and passing through (x, y). Return the discriminant and the solutions (z1, z2).
     * 
     * If there are no solutions, the solutions (z1, z2) are still computed as if the discriminant
     * was zero. In general you should ignore these values if the discriminant is negative. However,
     * if you know that the discriminant is ideally nonnegative but floating-point imprecision may
     * produce a negative discriminant, these values are needed. I found that this produces simpler
     * code than using an Option for when there are no solutions.
     */
    pub fn intersect_z_line_core(&self, xy: &Point2<f64>) -> (f64, f64, f64) {
        let m = self.matrix();
        let x = xy.x;
        let y = xy.y;
        let a1 = m[(0, 2)]; 
        let b1 = m[(0, 0)] * x + m[(0, 1)] * y + m[(0, 3)]; 
        let a2 = m[(1, 2)]; 
        let b2 = m[(1, 0)] * x + m[(1, 1)] * y + m[(1, 3)]; 
        let a = squared(a1) + squared(a2);
        let b = 2.0 * (a1 * b1 + a2 * b2);
        let c = squared(b1) + squared(b2) - 1.0;
        let discriminant = squared(b) - 4.0 * a * c;
        let discriminant_clipped = discriminant.max(0.0);
        let tmp = 1.0 / (2.0 * a);
        let z1 = (-b - discriminant_clipped.sqrt()) * tmp;
        let z2 = (-b + discriminant_clipped.sqrt()) * tmp;
        (discriminant, z1, z2)
    }

    /**
     * Return true if the line (cos(theta), sin(theta), z) intersects the pipe section.
     */
    pub fn intersects_z_line_theta(&self, theta: f64) -> bool {
        let p = Point2::new(theta.cos(), theta.sin());
        self.intersect_z_line_core(&p).0 >= 0.0
    }

    /**
     * Intersect the line (cos(theta), sin(theta), z). It is assumed that there is an intersection.
     */
    pub fn intersect_z_line_theta(&self, theta: f64) -> (Point3<f64>, Point3<f64>) {
        let p = Point2::new(theta.cos(), theta.sin());
        let (_, z1, z2) = self.intersect_z_line_core(&p);
        (
            Point3::new(p.x, p.y, z1),
            Point3::new(p.x, p.y, z2),
        )
    }

    /**
     * Project the cylinder onto the plane z = 0, producing a stripe whose boundary is two
     * parallel lines. Return these lines.
     */
    fn z_plane_projection_boundary(&self) -> (Line2D, Line2D) {
        let (p_3d, d_3d) = self.axis_line();
        // Intersection of axial line with the z = 0 plane.
        let p = Point2::new(p_3d.x, p_3d.y);
        // Projection of the direction of the axial line onto the z = 0 plane.
        let d = Vector2::new(d_3d.x, d_3d.y).normalize();
        // Unit vector orthogonal to d.
        let d_ortho = Vector2::new(-d.y, d.x);
        // 2x2 rotation matrix that orients the stripe so it is parallel to the x-axis.
        // If this matrix is R, Rd = [1, 0].
        let rotation = Matrix2::new(
            d.x, d.y,
            -d.y, d.x,
        );
        let corrected_ellipse_matrix = rotation * self.inv_top_left_matrix();
        let half_stripe_width = f64::hypot(
            corrected_ellipse_matrix[(1, 0)],
            corrected_ellipse_matrix[(1, 1)]
        );
        (
            Line2D { p: p - d_ortho * half_stripe_width, d },
            Line2D { p: p + d_ortho * half_stripe_width, d }
        )
    }

    /**
     * Let K be the intersection of this Cylinder with the "base cylinder" x^2 + y^2 = 1. The
     * projection of K into the plane z = 0 is a subset of the circle (x, y, 0) with x^2 + y^2 = 1.
     * There are four possibilities with this projection:
     * 
     * 1. It is empty.
     * 2. It is the entire circle.
     * 3. It is a single closed arc in the circle.
     * 4. It is two non-overlapping closed arcs in the circle.
     * 
     * The 2 or 4 bounding values of the closed arc or arcs are called the "critical thetas," hence
     * the name of this method.
     * 
     * The four possible cases are described by the CylinderIntersectionSolutions enum returned by this
     * method.
     */
    pub fn get_critical_thetas(&self) -> CylinderIntersectionSolutions {
        let (line1, line2) = self.z_plane_projection_boundary();
        let mut solutions1 = line1.intersect_unit_circle();
        let mut solutions2 = line2.intersect_unit_circle();

        if let (Some((p1, p2)), Some((p3, p4))) = (solutions1, solutions2) {
            let (t1, t2, t3, t4) = sort4((angle(&p1), angle(&p2), angle(&p3), angle(&p4)));

            // Special case: if t1 == t2, then (t1 + t2) / 2 will always be in K. Instead test the
            // [t2, t3] interval and invert the result.
            let tolerance = 1e-5;
            let (t_test, flip) = if (t1 - t2).abs() < tolerance {
                ((t2 + t3) / 2.0, true)
            } else {
                ((t1 + t2) / 2.0, false)
            };

            if self.intersects_z_line_theta(t_test) != flip {
                return CylinderIntersectionSolutions::TwoIntervals(
                    CircularInterval::new(t1, t2),
                    CircularInterval::new(t3, t4),
                );
            } else {
                return CylinderIntersectionSolutions::TwoIntervals(
                    CircularInterval::new(t2, t3),
                    CircularInterval::new(t4, t1 + f64::consts::TAU),
                );
            }
        }

        if let (None, Some(_)) = (solutions1, solutions2) {
            (solutions1, solutions2) = (solutions2, solutions1);
        }

        if let (Some((p1, p2)), None) = (solutions1, solutions2) {
            let (t1, t2) = sort2((angle(&p1), angle(&p2)));
            let mid_angle = (t1 + t2) / 2.0;
            if self.intersects_z_line_theta(mid_angle) {
                return CylinderIntersectionSolutions::OneInterval(
                    CircularInterval::new(t1, t2)
                );
            } else {
                return CylinderIntersectionSolutions::OneInterval(
                    CircularInterval::new(t2, t1 + f64::consts::TAU)
                );
            }
        }

        if self.intersects_z_line_theta(0.0) {
            CylinderIntersectionSolutions::All
        } else {
            CylinderIntersectionSolutions::Empty
        }
    }

    /// Intersect this cylinder with a plane at z, parallel to the xy-plane. Return the point on the
    /// resulting ellipse parametrized by angle theta from 0 to 2pi.
    pub fn intersect_z_plane_parametrized(&self, z: f64, theta: f64) -> Point3<f64> {
        let ellipse_center = self.intersect_axis_line_z_plane(z);
        let xy = Vector2::new(theta.cos(), theta.sin());
        let displacement_2d = self.inv_top_left_matrix() * xy;
        ellipse_center + Vector3::new(displacement_2d.x, displacement_2d.y, 0.0)
    }

    /// Inverse of intersect_z_plane_parametrized, returning the "theta" value.
    pub fn z_plane_ellipse_theta(&self, z: f64, p: &Point2<f64>) -> f64 {
        let ellipse_center = self.intersect_axis_line_z_plane(z).xy();
        let displacement_2d = p - ellipse_center;
        let xy = self.top_left_matrix() * displacement_2d;
        angle_vector(&xy)
    }

    /// Intersect this cylinder with a plane at z, parallel to the xy-plane. Return the result as a
    /// CCurve.
    pub fn intersect_z_plane(&self, z: f64) -> CCurve {
        CCurve::plane(self.clone(), z)
    }

    /// Intersect this cylinder with a planes at +z and -z, parallel to the xy-plane. Return the
    /// result as two CCurves.
    pub fn intersect_z_planes(&self, z: f64) -> CylinderIntersection {
        CylinderIntersection {
            ccurves: vec![
                self.intersect_z_plane(z),
                self.intersect_z_plane(-z),
            ]
        }
    }

    /**
     * Intersect this Cylinder with the base Cylinder and return the connected components as a set
     * of CCurves.
     */
    pub fn intersect_base_cylinder(&self) -> CylinderIntersection {
        let solutions = self.get_critical_thetas();
        let ccurves = match solutions {
            CylinderIntersectionSolutions::Empty => vec![],
            CylinderIntersectionSolutions::All => {
                vec![
                    CCurve::wrapped_loop(self.clone(), false),
                    CCurve::wrapped_loop(self.clone(), true),
                ]
            },
            CylinderIntersectionSolutions::OneInterval(interval) => {
                vec![
                    CCurve::side_loop(self.clone(), interval.start, interval.end),
                ]
            },
            CylinderIntersectionSolutions::TwoIntervals(interval1, interval2) => {
                vec![
                    CCurve::side_loop(self.clone(), interval1.start, interval1.end),
                    CCurve::side_loop(self.clone(), interval2.start, interval2.end),
                ]
            },
        };
        CylinderIntersection { ccurves }
    }

    /**
     * Intersect this Cylinder with another Cylinder.
     */
    pub fn intersect_cylinder(&self, other: &Cylinder) -> CylinderIntersection {
        // Let D(M_1) be self and let D(M_2) be other.
        // Note that D(M) = M^-1 D(I), so:
        //
        // intersect(D(M_1), D(M_2)) = M_2^-1 intersect(D(M_1 M_2^-1), D(I))
        //
        // Below, transformed_cylinder is D(M_1 M_2^-1) which we intersect with
        // the base cylinder, and M_2^-1 is applied to the result.
        let transform = other.transformation_from_base_cylinder();
        let transformed_cylinder = Cylinder::from_matrix_unchecked(
            self.matrix() * other.inv_matrix()
        );

        let untransformed_intersection = transformed_cylinder.intersect_base_cylinder();

        untransformed_intersection.transform(&transform)
    }
}


#[cfg(test)]
mod test {
    use core::f64;
    use std::io::pipe;

    use approx::*;
    use crate::mesh::{ColoredMesh, MeshLike, Color};
    use crate::{cylinder, mesh::Mesh};
    use crate::cylinder_curve::CCurveKind;
    use crate::config::{CylinderMeshConfig, TorusMeshConfig};
    use std::path::PathBuf;

    use super::*;
    use na::{Vector2, Point3};

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

    fn example_cylinder_2() -> Cylinder {
        Cylinder {
            m11: -1.3,
            m12: -0.1,
            m13: 0.3,
            m14: -0.6,
            m21: 0.3,
            m22: -1.4,
            m23: 0.6,
            m24: -0.1,
        }
    }

    #[test]    
    fn test_intersect_line_circle_1() {
        let line = Line2D {
            p: Point2::new(0.0, 0.0),
            d: Vector2::new(1.0, 1.0).normalize(),
        };
        if let Some(points) = line.intersect_unit_circle() {
            let tmp = 0.5f64.sqrt();
            assert_abs_diff_eq!(points.0, Point2::new(tmp, tmp));
            assert_abs_diff_eq!(points.1, Point2::new(-tmp, -tmp));
        } else {
            panic!("Line doesn't intersect circle");
        }
    }

    #[test]    
    fn test_intersect_line_circle_no_intersection() {
        let line = Line2D {
            p: Point2::new(2.0, 0.0),
            d: Vector2::new(1.0, 1.0).normalize(),
        };
        assert!(matches!(line.intersect_unit_circle(), None));
    }

    #[test]
    fn test_intersect_z_line() {
        let cylinder = example_cylinder();
        let x = 0.0;
        let y = 0.0;
        let xy = Point2::new(x, y);
        let (d, z1, z2) = cylinder.intersect_z_line_core(&xy);
        assert!(d > 0.0);
        assert_abs_diff_eq!(cylinder.scalar_field(&Point3::new(x, y, z1)), 0.0);
        assert_abs_diff_eq!(cylinder.scalar_field(&Point3::new(x, y, z2)), 0.0);
    }

    /**
     * Lines erected from the z-plane projection boundary should have a discriminant of about 0.
     */
    #[test]
    fn test_z_plane_projection_boundary() {
        let cylinder = example_cylinder();
        let (line1, line2) = cylinder.z_plane_projection_boundary();
        for line in [line1, line2] {
            for t in [0.0, 1.2, -3.0] {
                let p = line.p + line.d * t;
                assert_abs_diff_eq!(
                    cylinder.intersect_z_line_core(&p).0,
                    0.0,
                    epsilon = 1e-10
                );
            }
        }
    }

    #[test]
    fn test_critical_thetas() {
        let cylinder = example_cylinder();
        let thetas = cylinder.get_critical_thetas().values();
        for theta in thetas.iter() {
            let p = Point2::new(theta.cos(), theta.sin());
            let (discriminant, _, _) = cylinder.intersect_z_line_core(&p);
            assert_abs_diff_eq!(discriminant, 0.0, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_strip_curve() {
        let cylinder = example_cylinder();
        let solutions = cylinder.get_critical_thetas();
        let n = 30;
        for i in 0..n {
            let theta = (i as f64) * f64::consts::TAU / (n as f64);
            let p = Point2::new(theta.cos(), theta.sin());
            let (discriminant, _, _) = cylinder.intersect_z_line_core(&p);
            let contains = solutions.contains(theta);
            if contains {
                assert!(discriminant >= 0.0);
            } else {
                assert!(discriminant < 0.0);
            }
        }
    }

    #[test]
    fn test_intersect_base_cylinder() {
        let cylinder = example_cylinder();
        let intersection = cylinder.intersect_base_cylinder();
        for curve in intersection.ccurves {
            let point = curve.at(0.25);
            assert_abs_diff_eq!(cylinder.scalar_field(&point), 0.0, epsilon = 1e-10);
            assert_abs_diff_eq!(point.x.hypot(point.y), 1.0, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_intersect_cylinder() {
        let cylinder = example_cylinder();
        let cylinder2 = example_cylinder_2();
        let intersection = cylinder.intersect_cylinder(&cylinder2);
        for curve in intersection.ccurves {
            let n = 30;
            for i in 0..n {
                let t = i as f64 / n as f64;
                let point = curve.at(t);
                assert_abs_diff_eq!(cylinder.scalar_field(&point), 0.0, epsilon = 1e-10);
                assert_abs_diff_eq!(cylinder2.scalar_field(&point), 0.0, epsilon = 1e-10);
            }
        }
    }

    #[test]
    fn test_intersect_z_plane() {
        let cylinder = example_cylinder();
        let z = 0.3;
        let curve = cylinder.intersect_z_plane(0.3);
        for t in [0.0, 0.5, 0.7] {
            let point = curve.at(t);
            assert_abs_diff_eq!(cylinder.scalar_field(&point), 0.0, epsilon = 1e-10);
            assert_abs_diff_eq!(point.z, z, epsilon = 1e-10);
        }
    }

    /// Bug found in cylinder intersections.
    #[test]
    fn test_intersect_bug() {
        let cylinder_1 = Cylinder {
            m11: 0.5,
            m12: -0.866,
            m13: 0.517,
            m14: 0.0,
            m21: -0.866,
            m22: -0.5,
            m23: 0.0,
            m24: -0.05,
        };
        let cylinder_2 = Cylinder {
            m11: 0.517,
            m12: 0.0,
            m13: 1.0,
            m14: 0.0,
            m21: 0.0,
            m22: -0.517,
            m23: 0.0,
            m24: -0.1,
        };
        let intersection = cylinder_1.intersect_cylinder(&cylinder_2);
        for curve in intersection.ccurves {
            for t in linspace(0.0, 1.0, 500) {
                let p = curve.at(t);
                assert_abs_diff_eq!(cylinder_1.scalar_field(&p), 0.0, epsilon = 1e-5);
                assert_abs_diff_eq!(cylinder_2.scalar_field(&p), 0.0, epsilon = 1e-5);
            }
        }
    }

    #[test]
    fn test_intersect_bug_2() {
        let cylinder = Cylinder {
            m11: 0.96,
            m12: 1.67,
            m13: -0.45,
            m14: 0.16,
            m21: -1.67,
            m22: 0.967,
            m23: 1.67,
            m24: 0.046,
        };
        let intersection = cylinder.intersect_base_cylinder();
        for curve in intersection.ccurves {
            assert!(matches!(curve.kind, CCurveKind::SideLoop(_, _)));
            for t in linspace(0.0, 1.0, 500) {
                let p = curve.at(t);
                assert_abs_diff_eq!(cylinder.scalar_field(&p), 0.0, epsilon = 1e-5);
            }
        }
    }

    #[test]
    #[ignore]
    fn test_intersect_bug_3() {
        let pipe_section_1 = PipeSection {
            a: 0.5000000000000001,
            b: -0.8660254037844386,
            c: 0.5176380902050415,
            d: 0.0,
            w: 0.0,
        };
        let pipe_section_2 = PipeSection {
            a: 0.5000000000000001,
            b: 0.8660254037844386,
            c: 0.5176380902050415,
            d: 0.0,
            w: 0.0,
        };
        let orthogonal_pipe_section = PipeSection {
            a: -0.5000000000000003,
            b: 0.8660254037844392,
            c: 1.931851652578138,
            d: 0.0,
            w: 0.0,
        };
        let intersection = pipe_section_1.intersect(&pipe_section_2);
        let cylinder_config = CylinderMeshConfig {
            half_length: 5.0,
            resolution: 0.01,
        };
        let curve_config = TorusMeshConfig {
            radius: 0.05,
            linear_segments: 50,
            radial_segments: 10,
        };
        ColoredMesh {
            meshes: vec![
                (pipe_section_1.as_mesh(&cylinder_config), Color { red: 255, green: 255, blue: 255 }),
                (pipe_section_2.as_mesh(&cylinder_config), Color { red: 255, green: 255, blue: 255 }),
                (orthogonal_pipe_section.as_mesh(&cylinder_config), Color { red: 128, green: 128, blue: 255 }),
                (intersection.as_mesh(&curve_config), Color { red: 255, green: 255, blue: 255 }),
            ]
        } .write_ply_file(&PathBuf::from("cylinders.ply"));
    }
}