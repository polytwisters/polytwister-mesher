use crate::pipe_section::PipeSection;
use na::{Point2, Point3, Vector2};
use crate::utils::squared;

/**
 * A 2D line given by {p + dt | t in R} where p and d are in R^2 and d is a unit vector.
 */
pub struct Line2D {
    pub p: Vector2<f64>,
    pub d: Vector2<f64>,
}

impl Line2D {
    /**
     * Return the intersection of this line with the circle x^2 + y^2 = 0. Returns either a tuple of
     * two intersection points (the same point twice if the line is tangent), or None if the line
     * does not intersect the circle.
     */
    pub fn intersect_unit_circle(&self) -> Option<(Vector2<f64>, Vector2<f64>)> {
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
        Some((
            self.p + t1 * self.d,
            self.p + t2 * self.d,
        ))
    }
}

impl PipeSection {

    /**
     * Given a 2D point (x, y), intersect the pipe section with the line parallel to the z-axis
     * and passing through (x, y). Return None if there are no solutions, otherwise return a tuple
     * of two z values satisfying the equations.
     */
    fn intersect_z_line(&self, xy: &Point2<f64>) -> Option<(f64, f64)> {
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
        if discriminant < 0.0 {
            None
        } else {
            let tmp = 1.0 / (2.0 * a);
            let z1 = (-b - discriminant.sqrt()) * tmp;
            let z2 = (-b + discriminant.sqrt()) * tmp;
            Some((z1, z2))
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use na::{Vector2};

    #[test]    
    fn test_intersect_line_circle_1() {
        let line = Line2D {
            p: Vector2::new(0.0, 0.0),
            d: Vector2::new(1.0, 1.0).normalize(),
        };
        if let Some(points) = line.intersect_unit_circle() {
            let tmp = 0.5f64.sqrt();
            assert_abs_diff_eq!(points.0, Vector2::new(tmp, tmp));
            assert_abs_diff_eq!(points.1, Vector2::new(-tmp, -tmp));
        } else {
            panic!("Line doesn't intersect circle");
        }
    }

    #[test]    
    fn test_intersect_line_circle_no_intersection() {
        let line = Line2D {
            p: Vector2::new(2.0, 0.0),
            d: Vector2::new(1.0, 1.0).normalize(),
        };
        assert!(matches!(line.intersect_unit_circle(), None));
    }

    #[test]
    fn test_intersect_z_line() {
        let pipe = PipeSection { a: 1.2, b: -0.5, c: 0.4, d: 0.1, w: 0.1 };
        let x = 0.0;
        let y = 0.0;
        let xy = Point2::new(x, y);
        if let Some((z1, z2)) = pipe.intersect_z_line(&xy) {
            assert_abs_diff_eq!(pipe.scalar_field(&Point3::new(x, y, z1)), 0.0);
            assert_abs_diff_eq!(pipe.scalar_field(&Point3::new(x, y, z2)), 0.0);
        } else {
            panic!("No intersection");
        }
    }
}