use na::{Vector2};

/**
 * A 2D line given by {p + dt | t in R} where p and d are in R^2.
 */
pub struct Line2D {
    pub p: Vector2<f64>,
    pub d: Vector2<f64>,
}

impl Line2D {
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

mod test {
    use crate::strip_curves::*;
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

}