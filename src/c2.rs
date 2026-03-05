use std::ops::{Mul, Div};

use na::{Complex, Vector2, Vector4, ComplexField, Matrix2};

/// A complex 2-vector, i.e., a member of the vector space C^2.
#[derive(Clone, Copy, Debug)]
pub struct C2 {
    pub vec: Vector2<Complex<f64>>
}

impl C2 {
    /// Create the C^2 vector (x, y) from x and y.
    pub fn new(x: Complex<f64>, y: Complex<f64>) -> Self {
        Self {
            vec: Vector2::new(x, y)
        }
    }

    /// Create the C^2 vector (a + bi, c + di) for real a, b, c, d.
    pub fn from_components(a: f64, b: f64, c: f64, d: f64) -> Self {
        Self {
            vec: Vector2::new(Complex::new(a, b), Complex::new(c, d))
        }
    }

    /// Create the C^2 vector (a + bi, c + di) from the R^4 vector (a, b, c, d).
    pub fn from_vector4(vec: &Vector4<f64>) -> Self {
        Self {
            vec: Vector2::new(
                Complex::new(vec.x, vec.y),
                Complex::new(vec.z, vec.w),
            )
        }
    }

    /// Convert the C^2 vector (a + bi, c + di) to the R^4 vector (a, b, c, d).
    pub fn to_vector4(&self) -> Vector4<f64> {
        Vector4::new(
            self.vec.x.re,
            self.vec.x.im,
            self.vec.y.re,
            self.vec.y.im,
        )
    }

    pub fn abs(&self) -> f64 {
        self.vec.norm()
    }

    pub fn abs_difference(&self, other: &Self) -> f64 {
        (self.vec - other.vec).norm()
    }

    pub fn inner(&self, other: &Self) -> Complex<f64> {
        other.vec.dotc(&self.vec)
    }

    pub fn inner_abs(&self, other: &Self) -> f64 {
        self.inner(&other).abs()
    }

    pub fn rotate_real_b(&self) -> Self {
        let tmp = self.vec.y.conj() / self.vec.y.abs();
        C2 { vec: self.vec * tmp }
    }

    /// Return the special 2x2 unitary matrix that sends this C2 to k[1, 0] for real k.
    fn normalizing_su2_matrix(&self) -> Matrix2<Complex<f64>> {
        let norm = self.abs();
        Matrix2::new(
            self.vec.x.conj() / norm,
            self.vec.y.conj() / norm,
            -self.vec.y / norm,
            self.vec.x / norm,
        )
    }

    /// Return the inverse of self.normalizing_su2_matrix.
    fn normalizing_su2_matrix_inv(&self) -> Matrix2<Complex<f64>> {
        self.normalizing_su2_matrix().adjoint()
    }

    fn intersect_pipes_core(pipe1: &Self, pipe2: &Self) -> Option<(C2, C2)> {
        let a1 = pipe1.vec.x;
        let b1 = pipe1.vec.y;
        let a2 = pipe2.vec.x;
        let b2 = pipe2.vec.y;

        let r1 = 1.0 / b1.re;
        let r2 = 1.0 / b2.re;
        let c1 = -a1.conj() / b1;
        let c2 = -a2.conj() / b2;
        let d = (c2 - c1).abs();
        let ell = (r1 * r1 - r2 * r2 + d * d) / (2.0 * d);
        let discriminant = r1 * r1 - ell * ell;
        if discriminant < 0.0 {
            return None;
        }
        let tmp = Complex::new(ell, -discriminant.sqrt());
        let z2_a = c1 + tmp / d * (c2 - c1);
        let z2_b = c1 + tmp.conj() / d * (c2 - c1);
        let solution_a = C2::new(Complex::new(1.0, 0.0), z2_a);
        let solution_b = C2::new(Complex::new(1.0, 0.0), z2_b);

        Some((solution_a, solution_b))
    }

    pub fn intersect_pipes(pipe1: &Self, pipe2: &Self, pipe3: &Self) -> Option<(C2, C2)> {
        let u = pipe3.normalizing_su2_matrix();
        let u_inv = pipe3.normalizing_su2_matrix_inv();
        let k = 1.0 / pipe3.abs();
        let pipe1_transformed = C2 { vec: (u * pipe1.vec) * Complex::from_real(k) } .rotate_real_b();
        let pipe2_transformed = C2 { vec: (u * pipe2.vec) * Complex::from_real(k) } .rotate_real_b();
        match Self::intersect_pipes_core(&pipe1_transformed, &pipe2_transformed) {
            Some((a, b)) => {
                Some((
                    C2 { vec: u_inv * a.vec * Complex::from_real(k) },
                    C2 { vec: u_inv * b.vec * Complex::from_real(k) }
                ))
            },
            None => None
        }
    }

    pub fn similarity(&self, other: &Self) -> f64 {
        self.inner_abs(&other) / (self.abs() * other.abs())
    }
}

impl Mul<f64> for &C2 {
    type Output = C2;
    fn mul(self, rhs: f64) -> C2 {
        C2 {
            vec: self.vec * Complex::from_real(rhs)
        }
    }
}

impl Div<f64> for &C2 {
    type Output = C2;
    fn div(self, rhs: f64) -> C2 {
        C2 {
            vec: self.vec / Complex::from_real(rhs)
        }
    }
}

#[cfg(test)]
mod test {
    use approx::*;
    use super::*;

    #[test]
    fn test_from_components() {
        let (a, b, c, d) = (1.0, 2.0, -3.0, 4.0);
        let x = C2::from_components(a, b, c, d);
        assert_eq!(
            x.vec,
            Vector2::new(
                Complex { re: a, im: b },
                Complex { re: c, im: d },
            )
        );
    }

    #[test]
    fn test_from_vector4() {
        let (a, b, c, d) = (1.0, 2.0, -3.0, 4.0);
        let vector4 = Vector4::<f64>::new(a, b, c, d);
        let x = C2::from_vector4(&vector4);
        assert_eq!(
            x.vec,
            Vector2::new(
                Complex { re: a, im: b },
                Complex { re: c, im: d },
            )
        );
    }

    #[test]
    fn test_rotate_real_b() {
        let x = C2::from_components(0.0, 0.3, 1.0, 0.3);
        let y = x.rotate_real_b();
        assert_abs_diff_eq!(y.vec.y.im, 0.0);
        assert_abs_diff_eq!(y.abs(), x.abs());
    }

    #[test]
    fn test_normalizing_su2_matrix() {
        let x = C2::from_components(0.1, 0.3, 1.0, 0.3);
        let m = x.normalizing_su2_matrix();
        let actual = m * x.vec;
        let expected = Vector2::new(Complex::new(x.abs(), 0.0), Complex::ZERO);
        assert_abs_diff_eq!((actual - expected).norm(), 0.0);
    }

    #[test]
    fn test_normalizing_su2_matrix_is_unitary() {
        let x = C2::from_components(0.1, 0.3, 1.0, 0.3);
        let m = x.normalizing_su2_matrix();
        let adjoint = m.adjoint();
        let actual = m * adjoint;
        let expected = Matrix2::identity();
        assert_abs_diff_eq!((actual - expected).norm(), 0.0);
    }

    #[test]
    fn test_intersect_3_pipes_core() {
        let pipe1 = C2::from_components(0.0, 0.3, 1.0, 0.3).rotate_real_b();
        let pipe2 = C2::from_components(1.0, 0.3, 1.0, 0.1).rotate_real_b();
        let intersection = C2::intersect_pipes_core(&pipe1, &pipe2);
        if let Some((p1, p2)) = intersection {
            for pipe in [pipe1, pipe2] {
                assert_abs_diff_eq!(p1.inner_abs(&pipe), 1.0);
                assert_abs_diff_eq!(p2.inner_abs(&pipe), 1.0);
            }
        } else {
            panic!("Didn't intersect");
        }
    }

    #[test]
    fn test_intersect_3_pipes() {
        let pipe1 = C2::from_components(0.0, 0.3, 1.0, 0.3);
        let pipe2 = C2::from_components(1.0, 0.3, 1.0, 0.1);
        let pipe3 = C2::from_components(1.1, -0.3, 0.4, 0.5); 
        let intersection = C2::intersect_pipes(&pipe1, &pipe2, &pipe3);
        if let Some((p1, p2)) = intersection {
            for pipe in [pipe1, pipe2, pipe3] {
                assert_abs_diff_eq!(p1.inner_abs(&pipe), 1.0, epsilon = 1e-5);
                assert_abs_diff_eq!(p2.inner_abs(&pipe), 1.0, epsilon = 1e-5);
            }
        } else {
            panic!("Didn't intersect");
        }
    }
}