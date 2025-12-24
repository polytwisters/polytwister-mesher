extern crate nalgebra as na;
use na::{Affine3, Matrix2, Matrix3, Matrix4, Rotation3, Vector2, Vector3, Point3};
use crate::utils::{squared};
use crate::ellipse_spacing::{warp_elliptic_angle, ellipse_circumference};

/**
 * An infinite hollow cylinder in R^3 created as the invertible affine
 * transformation of the "base cylinder" { (x, y, z): x^2 + y^2 = 1 }.
 */
#[derive(Clone, Copy, Debug)]
pub struct Cylinder {
    pub m11: f64,
    pub m12: f64,
    pub m13: f64,
    pub m14: f64,
    pub m21: f64,
    pub m22: f64,
    pub m23: f64,
    pub m24: f64,
}


impl Cylinder {
    pub fn base() -> Self {
        Cylinder { m11: 1.0, m12: 0.0, m13: 0.0, m14: 0.0, m21: 0.0, m22: 1.0, m23: 0.0, m24: 0.0 }
    }

    pub fn example() -> Cylinder {
        Cylinder {
            m11: 0.1,
            m12: -0.5,
            m13: 0.4,
            m14: -0.1,
            m21: 0.5,
            m22: -1.2,
            m23: 0.2,
            m24: 0.05,
        }
    }

    /**
     * Given a 4x4 matrix, use its top two rows as the cylinder's matrix. It is not checked that the
     * bottom half is the same as the bottom half of an identity matrix.
     */
    pub fn from_matrix_unchecked(matrix: Matrix4<f64>) -> Self {
        Cylinder {
            m11: matrix.m11,
            m12: matrix.m12,
            m13: matrix.m13,
            m14: matrix.m14,
            m21: matrix.m21,
            m22: matrix.m22,
            m23: matrix.m23,
            m24: matrix.m24 
        }
    }

    /**
     * Evaluate the scalar field:
     * 
     *     F(x, y, z) = (m_11 x + m_12 y + m_13 z + m_14)^2 + (m_21 x + m_22 y + m_23 z + m_24)^2 - 1
     * 
     * The cylinder is given by the isosurface F(x, y, z) = 0. For other points,
     * if F(x, y, z) < 0 then the point is inside the cylinder, F(x, y, z) > 0
     * is outside the cylinder, and F(x, y, z) = -1 is on the cylinder's symmetry
     * axis (plane of symmetry if a = b = 0).
     */
    pub fn scalar_field(&self, point: &Point3<f64>) -> f64 {
        squared(self.m11 * point.x + self.m12 * point.y + self.m13 * point.z + self.m14)
        + squared(self.m21 * point.x + self.m22 * point.y + self.m23 * point.z + self.m24)
        - 1.0
    }

    /**
     * Return an explicit parametrization of the line which is this cylinder section's symmetry axis.
     * The parametrization is v(t) = v_0 + d * t and returned as (v_0, d) so that v_0 is the
     * starting point and d is the direction vector. d is always a unit vector.
     */
    pub fn axis_line(&self) -> (Point3<f64>, Vector3<f64>) {
        // Solve M(x, y, z, 1) = 0 with z = 0 and z = 1 respectively. Ignore bottom two rows and
        // (x, y) = -inv_top_left (z, 1)
        let tmp = self.inv_top_left_matrix();
        let tmp2 = self.top_right_matrix();
        let point1_2d = -tmp * (tmp2 * Vector2::new(0.0, 1.0));
        let point2_2d = -tmp * (tmp2 * Vector2::new(1.0, 1.0));
        let point1 = Point3::new(point1_2d.x, point1_2d.y, 0.0);
        let point2 = Point3::new(point2_2d.x, point2_2d.y, 1.0);
        let d = (point2 - point1).normalize();
        (point1, d)
    }

    /**
     * The 4x4 matrix M:
     * 
     * M = [
     *     [m_11 m_12 m_13 m_14]
     *     [m_21 m_22 m_23 m_24]
     *     [0    0    1    0   ]
     *     [0    0    0    1   ]
     * ]
     * 
     * which defines the cylinder as |M(x, y, z, 1)|^2 = 1.
     */
    pub fn matrix(&self) -> Matrix4<f64> {
        Matrix4::new(
            self.m11, self.m12, self.m13, self.m14,
            self.m21, self.m22, self.m23, self.m24,
            0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 1.0,
        )
    }

    /**
     * The inverse of the top left 2x2 entries of M, defined in Cylinder::matrix(). Reused in
     * several places.
     */
    pub fn inv_top_left_matrix(&self) -> Matrix2<f64> {
        let m = self.matrix();
        Matrix2::new(
            m[(0, 0)], m[(0, 1)],
            m[(1, 0)], m[(1, 1)],
        ).try_inverse().unwrap()
    }

    /**
     * Top right 2x2 entries of M, defined in Cylinder::matrix(). Reused in several places.
     */
    fn top_right_matrix(&self) -> Matrix2<f64> {
        let m = self.matrix();
        Matrix2::new(
            m[(0, 2)], m[(0, 3)],
            m[(1, 2)], m[(1, 3)],
        )
    }

    /**
     * Inverse of Cylinder::matrix.
     */
    pub fn inv_matrix(&self) -> Matrix4<f64> {
        // Split M into [[A B] [0 I]] where 0 is a 2x2 zero matrix and I is a 2x2 identity matrix.
        // Block matrix inversion gives M^-1 = [[A^-1 -A^-1 B] [0 I]].
        // Thus only a 2x2 matrix inversion is needed, fortunately.
        let tmp1 = self.inv_top_left_matrix();
        let tmp2 = -tmp1 * self.top_right_matrix();
        Matrix4::new(
            tmp1.m11, tmp1.m12, tmp2.m11, tmp2.m12,
            tmp1.m21, tmp1.m22, tmp2.m21, tmp2.m22,
            0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 1.0,
        )
    }
    
    /**
     * Define the base cylinder C = {(x, y, z): x^2 + y^2 = 1, z in R}. Return the affine
     * transformation T such that T(C) is this cylinder section.
     */
    pub fn transformation_from_base_cylinder(&self) -> Affine3<f64> {
        Affine3::from_matrix_unchecked(self.inv_matrix())
    }

    /**
     * Define the base cylinder C = {(x, y, z): x^2 + y^2 = 1, z in R}. Return the affine
     * transformation T such that applying T to this cylinder section produces C.
     */
    pub fn transformation_to_base_cylinder(&self) -> Affine3<f64> {
        Affine3::from_matrix_unchecked(self.matrix())
    }

    /**
     * Compute the two displacement vectors that form the elliptic cross-section of the cylinder.
     */
    fn ellipse_vertex_displacements(&self) -> (Vector3<f64>, Vector3<f64>) {
        let (_, direction) = self.axis_line();
        let z = Vector3::z();
        // R rotates the cylinder's axis of symmetry to (0, 0, 1).
        let r = Rotation3::rotation_between(&direction, &z).unwrap_or(Rotation3::identity());
        let r_inv = Rotation3::rotation_between(&z, &direction).unwrap_or(Rotation3::identity());
        let project_xy: Matrix3<f64> = Matrix3::from_diagonal(&Vector3::new(1.0, 1.0, 0.0));
        let x = Vector3::x();
        let y = Vector3::y();
        // We specifically use "transform_vector" and not "transform_point" and are deliberately
        // ignoring the translation part of the affine transformation.
        let da = r_inv * project_xy * r * self.transformation_from_base_cylinder().transform_vector(&x);
        let db = r_inv * project_xy * r * self.transformation_from_base_cylinder().transform_vector(&y);
        if db.norm_squared() > da.norm_squared() {
            (db, da)
        } else {
            (da, db)
        }
    }

    fn ellipse_vertices(&self) -> (Point3<f64>, Point3<f64>) {
        let (start, _) = self.axis_line();
        let (da, db) = self.ellipse_vertex_displacements();
        (start + da, start + db)
    }

    fn surface_coords_to_cartesian(&self, u: f64, theta: f64) -> Point3<f64> {
        let (start, direction) = self.axis_line();
        let (da, db) = self.ellipse_vertex_displacements();
        let a = da.norm();
        let b = db.norm();
        let theta_warped = warp_elliptic_angle(theta, a, b);
        let cx = theta_warped.cos();
        let cy = theta_warped.sin();
        start + u * direction + cx * da + cy * db
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use approx::*;

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
            m11: 1.3,
            m12: -0.4,
            m13: 0.1,
            m14: -0.5,
            m21: 0.2,
            m22: -1.5,
            m23: 0.1,
            m24: -0.3,
        }
    }

    #[test]
    fn test_cylinder_axis() {
        let cylinder = example_cylinder();
        let (start, direction) = cylinder.axis_line();
        assert_abs_diff_eq!(cylinder.scalar_field(&start), -1.0);
        let point_2 = start + direction;
        assert_abs_diff_eq!(cylinder.scalar_field(&point_2), -1.0);
        let point_3 = start + direction * 1.2345;
        assert_abs_diff_eq!(cylinder.scalar_field(&point_3), -1.0);
    }

    #[test]
    fn test_inv_matrix() {
        let cylinder = example_cylinder();
        assert_abs_diff_eq!(
            cylinder.matrix() * cylinder.inv_matrix(),
            Matrix4::identity()
        )
    }
}