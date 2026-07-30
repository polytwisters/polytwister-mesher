use na::Vector4;
use na::{Complex, Vector2, ComplexField, Matrix2};
use crate::c2::C2;
use crate::pipe_section::PipeSection;
use crate::ring::RingSection;

const EPSILON: f32 = 1e-5;

#[derive(Clone, Copy, Debug)]
pub struct Fiber {
    pub vec: C2
}

impl Fiber {
    pub fn new(vec: C2) -> Self { Self { vec } }
    pub fn from_vector4(vec4: &Vector4<f32>) -> Self { Self::new(C2::from_vector4(vec4)) }
    pub fn to_vector4(&self) -> Vector4<f32> { self.vec.to_vector4() }

    pub fn zero() -> Self { Self { vec: C2::zero() } }

    pub fn radius(&self) -> f32 {
        self.vec.abs()
    }

    pub fn scale(&self, k: f32) -> Self {
        Self::new(&self.vec * k)
    }

    pub fn similarity(&self, other: &Self) -> f32 {
        self.vec.similarity(&other.vec)
    }

    pub fn deduplicate(fibers: &Vec<Self>, epsilon: f32) -> Vec<Self> {
        let epsilon = 1e-5;
        let mut result = vec![];
        for fiber in fibers.iter() {
            if result.iter().all(|fiber2| fiber.similarity(fiber2) < 1.0 - epsilon) {
                result.push(fiber.clone());
            }
        }
        result
    }

    pub fn cross_section(&self, w: f32) -> RingSection {
        RingSection::from_vector4(&self.to_vector4(), w)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Pipe {
    pub vec: C2
}

impl Pipe {
    pub fn new(vec: C2) -> Self { Self { vec } }
    pub fn from_vector4(vec4: &Vector4<f32>) -> Self { Self::new(C2::from_vector4(vec4)) }
    pub fn to_vector4(&self) -> Vector4<f32> { self.vec.to_vector4() }

    pub fn cross_section(&self, w: f32) -> PipeSection {
        PipeSection::from_vector4(&self.to_vector4(), w)
    }

    pub fn scale(&self, k: f32) -> Self {
        Self::new(&self.vec / k)
    }

    pub fn inner_abs(&self, fiber: &Fiber) -> f32 {
        self.vec.inner_abs(&fiber.vec)
    }

    /// For pipe P(y), compute |<y, x>|^2 - 1. This is 0 on the pipe, negative inside, and positive
    /// outside.
    pub fn scalar_field(&self, point: &C2) -> f32 {
        self.vec.inner_abs_squared(point) - 1.0
    }

    pub fn contains(&self, point: &C2, epsilon: f32) -> bool {
        self.scalar_field(point).abs() < epsilon
    }

    /// Given p_1, p_2 in C^2, solve the system of nonlinear equations
    /// |<p_1, z>| = 1, |<p_2, z>| = 1 for all z of the form z = (1, z2). If there are no solutions
    /// return None, if there are two solutions return both, if there is one solution return that
    /// solution duplicated.
    fn intersect_core(pipe1: &C2, pipe2: &C2) -> Option<(C2, C2)> {
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

    /// Given three C^2 vectors p_1, p_2, p_3, solve for all values z such that |<z, p_i>| = 1 for
    /// i in (1, 2, 3), producing a system of three nonlinear equations. If there are no solutions,
    /// None is returned. Otherwise, two solutions z1 and z2 are returned such that all phase
    /// rotations of z1 and z2 constitute the solutions to the system of equations.
    /// 
    /// Using the P function defined in the polytwister paper, this computes the intersection of the
    /// pipes P(p_1), P(p_2), and P(p_3).
    pub fn intersect(pipe1: &Self, pipe2: &Self, pipe3: &Self) -> Option<(Fiber, Fiber)> {
        let vec1 = pipe1.vec;
        let vec2 = pipe2.vec;
        let vec3 = pipe3.vec;

        // Normalize the three pipes so that pipe3 becomes (1, 0), reducing the problem to
        // intersect_core.
        let u = vec3.normalizing_su2_matrix();
        let u_inv = vec3.normalizing_su2_matrix_inv();
        let k = 1.0 / vec3.abs();
        let pipe1_transformed = C2 { vec: (u * vec1.vec) * Complex::from_real(k) } .rotate_real_b();
        let pipe2_transformed = C2 { vec: (u * vec2.vec) * Complex::from_real(k) } .rotate_real_b();
        match Self::intersect_core(&pipe1_transformed, &pipe2_transformed) {
            Some((a, b)) => {
                // Undo the normalizing transformation.
                Some((
                    Fiber::new(C2 { vec: u_inv * a.vec * Complex::from_real(k) }),
                    Fiber::new(C2 { vec: u_inv * b.vec * Complex::from_real(k) }),
                ))
            },
            None => None
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Log {
    pub vec: C2
}

impl Log {
    pub fn new(vec: C2) -> Self { Self { vec } }
    pub fn from_vector4(vec4: &Vector4<f32>) -> Self { Self::new(C2::from_vector4(vec4)) }
    pub fn to_vector4(&self) -> Vector4<f32> { self.vec.to_vector4() }

    pub fn scale(&self, k: f32) -> Self {
        Self::new(&self.vec / k)
    }

    /// For log L(y), compute |<y, x>|^2 - 1. This is 0 on the boundary of the log, negative in its
    /// interior, and positive in the exterior of the log.
    pub fn scalar_field(&self, point: &C2) -> f32 {
        self.vec.inner_abs_squared(point) - 1.0
    }

    /// Return true if this log contains the given point.
    pub fn contains(&self, point: &C2, epsilon: f32) -> bool {
        self.scalar_field(&point) < epsilon
    }

    /** The pipe bounding this log. */
    pub fn pipe(&self) -> Pipe {
        Pipe::new(self.vec)
    }
}

#[cfg(test)]
mod test {
    use approx::*;
    use super::*;

    #[test]
    fn test_intersect_3_pipes_core() {
        let pipe1 = C2::from_components(0.0, 0.3, -1.5, 0.33).rotate_real_b();
        let pipe2 = C2::from_components(1.1, -0.4, 0.2, 0.1).rotate_real_b();
        let intersection = Pipe::intersect_core(&pipe1, &pipe2);
        if let Some((p1, p2)) = intersection {
            for pipe in [pipe1, pipe2] {
                assert_abs_diff_eq!(p1.inner_abs(&pipe), 1.0, epsilon = 1e-6);
                assert_abs_diff_eq!(p2.inner_abs(&pipe), 1.0, epsilon = 1e-6);
            }
        } else {
            panic!("Didn't intersect");
        }
    }

    #[test]
    fn test_intersect_3_pipes() {
        let pipe1 = Pipe::new(C2::from_components(0.0, 0.3, 1.0, 0.3));
        let pipe2 = Pipe::new(C2::from_components(-1.0, 0.4, 0.5, 0.1));
        let pipe3 = Pipe::new(C2::from_components(1.1, -0.3, 0.4, 0.5)); 
        let intersection = Pipe::intersect(&pipe1, &pipe2, &pipe3);
        if let Some((f1, f2)) = intersection {
            for pipe in [pipe1, pipe2, pipe3] {
                assert!(pipe.contains(&f1.vec, EPSILON));
                assert!(pipe.contains(&f2.vec, EPSILON));
            }
        } else {
            panic!("Didn't intersect");
        }
    }
}