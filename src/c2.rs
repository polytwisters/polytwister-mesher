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

    /// Convert (a + bi, c + di) to the 4-tuple of floats (a, b, c, d).
    pub fn to_components(&self) -> (f64, f64, f64, f64) {
        (
            self.vec.x.re,
            self.vec.x.im,
            self.vec.y.re,
            self.vec.y.im,
        )
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

    pub fn zero() -> Self {
        Self { vec: Vector2::zeros() }
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

    /// Return the norm: ||(x, y)|| = sqrt(|x|^2 + |y|^2).
    pub fn abs(&self) -> f64 {
        self.vec.norm()
    }

    /// Return the absolute difference between two C2 vectors.
    pub fn abs_difference(&self, other: &Self) -> f64 {
        (self.vec - other.vec).norm()
    }

    /// Take the inner product <(x1, x2), (y1, y2)> = x1 conj(y1) + x2 conj(y2), which is linear in
    /// the first argument.
    pub fn inner(&self, other: &Self) -> Complex<f64> {
        other.vec.dotc(&self.vec)
    }

    /// Absolute value of inner product squared.
    pub fn inner_abs_squared(&self, other: &Self) -> f64 {
        self.inner(&other).modulus_squared()
    }

    /// Absolute value of inner product.
    pub fn inner_abs(&self, other: &Self) -> f64 {
        self.inner(&other).abs()
    }

    /// Given this vector (x, y), return a new vector (kx, ky) where |k| = 1 such that y is a real
    /// number. This will produce NaNs if y = 0.
    pub fn rotate_real_b(&self) -> Self {
        let tmp = self.vec.y.conj() / self.vec.y.abs();
        C2 { vec: self.vec * tmp }
    }

    /// Given this vector x, return a special 2x2 unitary matrix M in SU(2) such that Mx = (k, 0)
    /// for real k.
    pub fn normalizing_su2_matrix(&self) -> Matrix2<Complex<f64>> {
        let norm = self.abs();
        Matrix2::new(
            self.vec.x.conj() / norm,
            self.vec.y.conj() / norm,
            -self.vec.y / norm,
            self.vec.x / norm,
        )
    }

    /// Return the inverse of self.normalizing_su2_matrix.
    pub fn normalizing_su2_matrix_inv(&self) -> Matrix2<Complex<f64>> {
        self.normalizing_su2_matrix().adjoint()
    }

    /// Given C^2 vectors z1 and z2, return the similarity |<z1, z2>| / (||z1|| ||z2||). This is a
    /// cosine-like similarity function which is 1 iff they are phase rotations of each other.
    pub fn similarity(&self, other: &Self) -> f64 {
        self.inner_abs(&other) / (self.abs() * other.abs())
    }
}


impl Mul<f64> for &C2 {
    type Output = C2;

    /// Multiply by a real scalar.
    fn mul(self, rhs: f64) -> C2 {
        C2 {
            vec: self.vec * Complex::from_real(rhs)
        }
    }
}

impl Div<f64> for &C2 {
    type Output = C2;

    /// Divide by a real scalar.
    fn div(self, rhs: f64) -> C2 {
        C2 {
            vec: self.vec / Complex::from_real(rhs)
        }
    }
}

#[cfg(test)]
mod test {
    use std::vec;

    use approx::*;
    use super::*;

    fn c2_example() -> C2 {
        C2::from_components(1.0, 0.2, -3.0, -0.1)
    }

    #[test]
    fn test_from_to_components() {
        let tuple4 = (1.0, 2.0, -3.0, 4.0);
        let (a, b, c, d) = tuple4;
        let x = C2::from_components(a, b, c, d);
        assert_eq!(
            x.vec,
            Vector2::new(
                Complex { re: a, im: b },
                Complex { re: c, im: d },
            )
        );
        assert_eq!(tuple4, x.to_components());
    }

    #[test]
    fn test_from_to_vector4() {
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
        assert_eq!(x.to_vector4(), vector4);
    }

    #[test]
    fn test_abs() {
        let x = C2::from_components(3.0, 4.0, 0.0, 0.0);
        assert_abs_diff_eq!(x.abs(), 5.0);
    }

    #[test]
    fn test_inner() {
        // <(i, 2i), (1, i)> = i conj(1) + 2i conj(i) = 2 + i
        // Test fails if "inner" is linear in the second argument instead of the first.
        let x = C2::from_components(0.0, 1.0, 0.0, 2.0);
        let y = C2::from_components(1.0, 0.0, 0.0, 1.0);
        let expected = Complex::new(2.0, 1.0);
        let inner = x.inner(&y);
        assert_eq!(inner, expected);
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
        assert_abs_diff_eq!((actual - expected).norm(), 0.0, epsilon = 1e-6);
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
}