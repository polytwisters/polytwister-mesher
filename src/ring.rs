use core::f64;

use na::{Complex, ComplexField};
use nalgebra::{Point2, Point3, Vector4};
use crate::polyline::{self, Polyline};
use crate::{mesh::Mesh, utils::angle};
use crate::config::{RingMeshConfig};

#[derive(Clone, Copy, Debug)]
pub struct RingSection {
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub d: f64,
    pub w: f64,
}

impl RingSection {
    pub fn from_vector4(vector: &Vector4<f64>, w: f64) -> Self {
        RingSection {
            a: vector.x,
            b: vector.y,
            c: vector.z,
            d: vector.w,
            w
        }
    }

    pub fn as_points(&self) -> RingSectionResult {
        let tolerance = 1e-10;
        let denom = self.c.hypot(self.d);
        if denom < tolerance {
            if self.w.abs() < tolerance {
                return RingSectionResult::XYCircle { radius: self.a.hypot(self.b) };
            } else {
                return RingSectionResult::Empty;
            }
        }
        let discriminant = self.w / denom;
        if discriminant > 1.0 {
            return RingSectionResult::Empty;
        }
        let theta = angle(&Point2::new(self.c, self.d));
        let phi1 = discriminant.asin();
        let phi2 = std::f64::consts::PI - phi1;
        let k1 = Complex::from_polar(1.0, phi1 - theta);
        let k2 = Complex::from_polar(1.0, phi2 - theta);
        let point1_c2 = (
            k1 * Complex::new(self.a, self.b),
            k1 * Complex::new(self.c, self.d),
        );
        let point2_c2 = (
            k2 * Complex::new(self.a, self.b),
            k2 * Complex::new(self.c, self.d),
        );
        RingSectionResult::Points(
            Point3::new(point1_c2.0.re, point1_c2.0.im, point1_c2.1.re),
            Point3::new(point2_c2.0.re, point2_c2.0.im, point2_c2.1.re),
        )
    }

    pub fn add_points_to_vec(&self, vec: &mut Vec<Point3<f64>>) {
        if let RingSectionResult::Points(p1, p2) = self.as_points() {
            vec.push(p1);
            vec.push(p2);
        }
    }

    pub fn as_mesh(&self, options: &RingMeshConfig) -> Mesh {
        match self.as_points() {
            RingSectionResult::Empty => {
                Mesh::empty()
            },
            RingSectionResult::Points(p1, p2) => {
                Mesh::merge(vec![
                    Mesh::uv_sphere(&p1, options.radius, options.longitudes, options.latitudes),
                    Mesh::uv_sphere(&p2, options.radius, options.longitudes, options.latitudes),
                ])
            },
            RingSectionResult::XYCircle { radius } => {
                let resolution = 500;
                let polyline = Polyline {
                    points: (0..resolution).map(|i| {
                        let t = i as f64 / resolution as f64;
                        let theta = t * f64::consts::TAU;
                        Point3::new(theta.cos(), theta.sin(), 0.0) * radius
                    }).collect::<_>(),
                    closed: true,
                };
                polyline.as_mesh(options.radius, options.longitudes)
            }
        }
    }
}

pub enum RingSectionResult {
    Empty,
    Points(Point3<f64>, Point3<f64>),
    XYCircle {
        radius: f64
    },
}

#[cfg(test)]
mod test {
    use super::*;
    use approx::*;

    #[test]
    fn test_ring_section() {
        let section = RingSection {
            a: 0.3,
            b: 0.4,
            c: -0.5,
            d: 0.1,
            w: 0.3,
        };
        assert!(matches!(section.as_points(), RingSectionResult::Points(_, _)));
    }

    #[test]
    fn test_ring_section_zero() {
        let section = RingSection {
            a: 1.0,
            b: 0.0,
            c: 0.0,
            d: 0.0,
            w: 0.0,
        };
        assert!(matches!(section.as_points(), RingSectionResult::XYCircle { radius: _ }));
    }
}