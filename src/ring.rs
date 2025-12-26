use na::{Complex, ComplexField};
use nalgebra::{Point2, Point3, Vector4};
use crate::{mesh::Mesh, utils::angle};

pub struct RingSection {
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub d: f64,
    pub w: f64,
}

pub struct RingMeshConfig {
    pub radius: f64,
    pub segments: usize,
    pub rings: usize,
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

    pub fn as_points(&self) -> Option<(Point3<f64>, Point3<f64>)> {
        let denom = self.c.hypot(self.d);
        if denom < 1e-10 {
            return None;
        }
        let discriminant = self.w / denom;
        if discriminant > 1.0 {
            return None;
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
        Some((
            Point3::new(point1_c2.0.re, point1_c2.0.im, point1_c2.1.re),
            Point3::new(point2_c2.0.re, point2_c2.0.im, point2_c2.1.re),
        ))
    }

    pub fn as_mesh(&self, options: &RingMeshConfig) -> Mesh {
        if let Some((center1, center2)) = self.as_points() {
            Mesh::merge(vec![
                Mesh::uv_sphere(&center1, options.radius, options.segments, options.rings),
                Mesh::uv_sphere(&center2, options.radius, options.segments, options.rings),
            ])
        } else {
            Mesh::empty()
        }
    }
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
        assert!(matches!(section.as_points(), Some(_)));
    }
}