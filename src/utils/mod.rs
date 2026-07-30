use core::f32;
use na::{Matrix2, Matrix4};
use nalgebra::{Point2, Vector2, Vector4};

mod unordered_triples;
pub use unordered_triples::UnorderedTriples;

pub mod adaptive_sampling;

pub fn squared(x: f32) -> f32 {
    x * x
}

pub fn angle(point: &Point2<f32>) -> f32 {
    f32::atan2(point.y, point.x).rem_euclid(f32::consts::TAU)
}

pub fn angle_vector(v: &Vector2<f32>) -> f32 {
    f32::atan2(v.y, v.x).rem_euclid(f32::consts::TAU)
}

pub fn sort2(x: (f32, f32)) -> (f32, f32) {
    let (x1, x2) = x;
    if x1 < x2 { (x1, x2) } else { (x2, x1) }
}

pub fn sort4(x: (f32, f32, f32, f32)) -> (f32, f32, f32, f32) {
    let mut tmp = [x.0, x.1, x.2, x.3];
    tmp.sort_by(f32::total_cmp);
    (tmp[0], tmp[1], tmp[2], tmp[3])
}

pub fn lerp(x1: f32, x2: f32, t: f32) -> f32 {
    x1 + (x2 - x1) * t
}

pub fn lerp_inverse(x1: f32, x2: f32, x: f32) -> f32 {
    (x - x1) / (x2 - x1)
}

/**
 * Return a linearly spaced series of n values in the closed interval [start, end]. Both endpoints
 * are inclusive.
 */
pub fn linspace(start: f32, end: f32, n: usize) -> Vec<f32> {
    (0..n).map(|i| start + (i as f32) / ((n - 1) as f32) * (end - start)).collect::<_>()
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
    pub matrix: Matrix2<f32>
}

impl Ellipse {
    /// Return 0.0 if the point p is on the ellipse, negative if inside, and positive if outside.
    fn scalar_field(&self, p: &Vector2<f32>) -> f32 {
        (self.matrix * p).norm_squared() - 1.0
    }

    /// Inverse of 2x2 matrix in the implicit equation. Its columns are vectors parallel to two
    /// conjugate diameters. They are not guaranteed orthogonal.
    fn inv_matrix(&self) -> Matrix2<f32> {
        self.matrix.try_inverse().unwrap()
    }

    /// Return two of the vertices of the ellipse, one along the major axis and one along the minor
    /// axis.
    /// 
    /// The vertices of an ellipse are the furthest and closest points from its center. Only two of
    /// the four vertices are returned. If the ellipse is a circle, the vertices are not defined,
    /// and this function returns any two points at 90-degree angles from each other.
    /// 
    /// The handedness convention is chosen so that the cross product v1 x v2 always points in the
    /// positive z direction.
    pub fn vertices(&self) -> (Vector2<f32>, Vector2<f32>) {
        let m = self.inv_matrix();
        let r1 = Vector2::new(m[(0, 0)], m[(1, 0)]);
        let r2 = Vector2::new(m[(0, 1)], m[(1, 1)]);
        // https://en.wikipedia.org/wiki/Ellipse#General_ellipse_2
        let denom = r1.norm_squared() - r2.norm_squared();
        let t = if denom.abs() <= f32::EPSILON {
            0.0
        } else {
            (2.0 * r1.dot(&r2) / denom).atan() / 2.0
        };
        let t2 = t + f32::consts::FRAC_PI_2;
        let v1 = r1 * t.cos() + r2 * t.sin();
        let v2 = r1 * t2.cos() + r2 * t2.sin();

        // Place the major axis first.
        let (major, minor) = if v1.norm_squared() > v2.norm_squared() {
            (v1, v2)
        } else {
            (v2, v1)
        };

        // To ensure correct handedness we compute the cross product. Its x- and y-coordinates are 0
        // because the two input vectors are orthogonal, and we only need its z-coordinate, which is
        // the determinant of the 2x2 matrix [major minor].
        let cross = major.x * minor.y - minor.x * major.y;
        let (major, minor) = if cross > 0.0 {
            (major, minor)
        } else {
            // To get a positive z-coordinate, we sign flip one of the vectors. Arbitrarily chose
            // the minor one.
            (major, -minor)
        };
        (major, minor)
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

/// Given a monotonic function f: [0, 1] -> bool, use bisection search to find an x in [0, 1] so
/// f(x) is right on the cusp of the switch from "true" to "false" or vice versa.
pub fn bisection_search<F: Fn(f32) -> bool>(f: F) -> f32 {
    let mut x_min = 0.0; 
    let mut x_max = 1.0;
    if f(x_min) == f(x_max) {
        return 1.0;
    }
    // True if the function is ramping from false to true.
    let upward = f(x_max);
    let mut x;
    for i in 0..5 {
        x = (x_min + x_max) / 2.0;
        // This flips the if/else statements if upward is false.
        if f(x) == upward {
            x_max = x;
        } else {
            x_min = x;
        }
    }
    (x_min + x_max) / 2.0
}

fn normalizing_su2_matrix(vec: &Vector4<f32>) -> Matrix4<f32> {
    let norm = vec.norm();
    // U = [
    //    [x* y*]
    //    [-y x]
    // ]
    // Use matrix representation: a + bi = [[a -b], [b a]]
    Matrix4::new(
        vec.x, vec.y, vec.z, vec.w,
        -vec.y, vec.x, -vec.w, vec.z,
        -vec.z, vec.w, vec.x, -vec.y,
        -vec.w, -vec.z, vec.y, vec.x,
    ) / norm
}

fn rotate_w_zero(vec: &Vector4<f32>) -> Vector4<f32> {
    let tmp = vec.z.hypot(vec.w);
    // Matrix: [[a -b], [b a]] [z, w] = [k, 0]
    let a = vec.z / tmp;
    let b = -vec.w / tmp;
    Vector4::new(
        a * vec.x - b * vec.y,
        b * vec.x + a * vec.y,
        a * vec.z - b * vec.w,
        b * vec.z + a * vec.w,
    )
}

pub fn torus_radius(pipe1: &Vector4<f32>, pipe2: &Vector4<f32>) -> f32 {
    let u = normalizing_su2_matrix(pipe2);
    let k = 1.0 / pipe2.norm();
    let p1n = rotate_w_zero(&(u * pipe1 * k));
    return 1.0f32.hypot((p1n.x.hypot(p1n.y) + 1.0) / p1n.z.abs()) * k;
}

#[cfg(test)]
mod test {
    use super::*;
    use approx::*;

    #[test]
    fn test_normalizing_su2_matrix() {
        let vec = Vector4::new(-2.3, 0.3, 1.4, -0.6);
        let u = normalizing_su2_matrix(&vec);
        assert_abs_diff_eq!(u * u.transpose(), Matrix4::identity());
        assert_abs_diff_eq!(u * vec, Vector4::new(vec.norm(), 0.0, 0.0, 0.0), epsilon = 1e-10);
    }

    #[test]
    fn test_rotate_real_w() {
        let vec = Vector4::new(-2.3, 0.3, 1.4, -0.6);
        let result = rotate_w_zero(&vec);
        assert_abs_diff_eq!(vec.norm(), result.norm());
        assert_abs_diff_eq!(result.w, 0.0);
    }
}