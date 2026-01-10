use core::f64;
use na::Matrix2;
use nalgebra::{Point2, Vector2};

pub fn squared(x: f64) -> f64 {
    x * x
}

pub fn angle(point: &Point2<f64>) -> f64 {
    f64::atan2(point.y, point.x).rem_euclid(f64::consts::TAU)
}

pub fn sort2(x: (f64, f64)) -> (f64, f64) {
    let (x1, x2) = x;
    if x1 < x2 { (x1, x2) } else { (x2, x1) }
}

pub fn sort4(x: (f64, f64, f64, f64)) -> (f64, f64, f64, f64) {
    let mut tmp = [x.0, x.1, x.2, x.3];
    tmp.sort_by(f64::total_cmp);
    (tmp[0], tmp[1], tmp[2], tmp[3])
}

pub fn lerp(x1: f64, x2: f64, t: f64) -> f64 {
    x1 + (x2 - x1) * t
}

/**
 * Return a linearly spaced series of n values in the closed interval [start, end]. Both endpoints
 * are inclusive.
 */
pub fn linspace(start: f64, end: f64, n: usize) -> Vec<f64> {
    (0..n).map(|i| start + (i as f64) / ((n - 1) as f64) * (end - start)).collect::<_>()
}

/**
 * Convert the vector of pairs: [(1, X), (2, 10), ..., (4, 6), (X, 5)] to the vector [1, 2, 3, ...].
 */
pub fn unzip_circle<T>(pairs: Vec<(T, T)>) -> Vec<T> {
    let (mut tmp1, mut tmp2): (Vec<T>, Vec<T>) = pairs.into_iter().unzip();
    tmp1.pop();
    tmp2.reverse();
    tmp2.pop();
    tmp1.append(&mut tmp2);
    tmp1
}

/// Ellipse in 2D space, centered on the origin, given by implicit equation:
/// (m11 x + m12 y)^2 + (m21 x + m22 y)^2 = 1.
pub struct Ellipse {
    pub matrix: Matrix2<f64>
}

impl Ellipse {
    /// Return 0.0 if the point p is on the ellipse, negative if inside, and positive if outside.
    fn scalar_field(&self, p: &Vector2<f64>) -> f64 {
        (self.matrix * p).norm_squared() - 1.0
    }

    /// Inverse of 2x2 matrix in the implicit equation. Its columns are vectors parallel to two
    /// conjugate diameters. They are not guaranteed orthogonal.
    fn inv_matrix(&self) -> Matrix2<f64> {
        self.matrix.try_inverse().unwrap()
    }

    /// Return two orthogonal vectors which are the parallel to the major and minor axes of the
    /// ellipse. If the ellipse is a circle and has no major or minor axes, any two orthogonal
    /// vectors on the circle are returned.
    pub fn vertices(&self) -> (Vector2<f64>, Vector2<f64>) {
        let m = self.inv_matrix();
        let r1 = Vector2::new(m[(0, 0)], m[(1, 0)]);
        let r2 = Vector2::new(m[(0, 1)], m[(1, 1)]);
        // Formula for the vertices of an ellipse.
        // https://en.wikipedia.org/wiki/Ellipse#General_ellipse_2
        let denom = r1.norm_squared() - r2.norm_squared();
        let t = if denom.abs() <= f64::EPSILON {
            0.0
        } else {
            (2.0 * r1.dot(&r2) / denom).atan() / 2.0
        };
        let t2 = t + f64::consts::FRAC_PI_2;
        let v1 = r1 * t.cos() + r2 * t.sin();
        let v2 = r1 * t2.cos() + r2 * t2.sin();
        if v1.norm_squared() > v2.norm_squared() {
            (v1, v2)
        } else {
            (v2, v1)
        }
    }
}

#[cfg(test)]
mod tests {
    use approx::*;
    use super::*;

    /// Construct an ellipse from known vertices. Ellipse::vertices() should return them.
    #[test]
    fn test_ellipse_vertices() {
        let v1 = Vector2::new(1.0, 1.3);
        let v2 = Vector2::new(-v1.y, v1.x) * 0.6;
        let inv_matrix = Matrix2::new(
            v1.x, v2.x,
            v1.y, v2.y,
        );
        let matrix = inv_matrix.try_inverse().unwrap();
        let ellipse = Ellipse { matrix };
        let (w1, w2) = ellipse.vertices();
        assert_abs_diff_eq!(v1, w1, epsilon = 1e-5);
        assert_abs_diff_eq!(v2, w2, epsilon = 1e-5);
    }

    #[test]
    fn test_ellipse_vertices_orthogonal() {
        let matrix = Matrix2::new(
            2.0, 3.0,
            6.0, 1.0,
        );
        let ellipse = Ellipse { matrix };
        let (v1, v2) = ellipse.vertices();
        assert_abs_diff_eq!(v1.dot(&v2), 0.0, epsilon = 1e-5);
        assert_abs_diff_eq!(ellipse.scalar_field(&v1), 0.0, epsilon = 1e-5);
        assert_abs_diff_eq!(ellipse.scalar_field(&v2), 0.0, epsilon = 1e-5);
    }
}