extern crate nalgebra as na;
use na::{Affine3, Matrix2, Matrix3, Matrix4, Rotation3, Vector2, Vector3};
use crate::mesh::{Mesh, Face};
use crate::utils::{squared};
use crate::ellipse_spacing::{warp_elliptic_angle, ellipse_circumference};

/**
 * A 3D cross section of a pipe. (PipeCrossSection felt too long.)
 */
#[derive(Clone, Copy, Debug)]
pub struct PipeSection {
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub d: f64,
    pub w: f64,
}


impl PipeSection {
    /**
     * Evaluate the scalar field associated with this pipe cross section:
     * 
     *     F(x, y, z) = (ax + by + cz + dw)^2 + (bx - ay + dz - cw)^2 - 1
     * 
     * The pipe cross section is given by the isosurface F(x, y, z) = 0. For other points, if
     * F(x, y, z) < 0 then the point is inside the pipe, F(x, y, z) > 0 is outside the pipe, and
     * F(x, y, z) = -1 is on the pipe's symmetry axis (plane of symmetry if a = b = 0).
     */
    pub fn scalar_field(&self, point: &Vector3<f64>) -> f64 {
        squared(self.a * point.x + self.b * point.y + self.c * point.z + self.d * self.w)
        + squared(self.b * point.x - self.a * point.y + self.d * point.z - self.c * self.w)
        - 1.0
    }

    /**
     * Return the pipe as a tuple (a, b, c, d, w).
     */
    fn abcdw(&self) -> (f64, f64, f64, f64, f64) {
        (self.a, self.b, self.c, self.d, self.w)
    }

    /**
     * Return an explicit parametrization of the line which is this pipe section's symmetry axis.
     * The parametrization is v(t) = v_0 + d * t and returned as (v_0, d) so that v_0 is the
     * starting point and d is the direction vector. d is always a unit vector.
     */
    fn axis_line(&self) -> (Vector3<f64>, Vector3<f64>) {
        // Solve M(x, y, z, 1) = 0 with z = 0 and z = 1 respectively. Ignore bottom two rows and
        // (x, y) = -inv_top_left (z, 1)
        let tmp = self.inv_top_left_matrix();
        let tmp2 = self.top_right_matrix();
        let point1_2d = -tmp * (tmp2 * Vector2::new(0.0, 1.0));
        let point2_2d = -tmp * (tmp2 * Vector2::new(1.0, 1.0));
        let point1 = Vector3::new(point1_2d.x, point1_2d.y, 0.0);
        let point2 = Vector3::new(point2_2d.x, point2_2d.y, 1.0);
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
     * which defines the pipe section as |M(x, y, z, w)|^2 = 1.
     */
    fn matrix(&self) -> Matrix4<f64> {
        let (a, b, c, d, w) = self.abcdw();
        Matrix4::new(
            a, b, c, d * w,
            b, -a, d, -c * w,
            0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 1.0,
        )
    }

    /**
     * The inverse of the top left 2x2 entries of M, defined in PipeSection::matrix(). Reused in
     * several places.
     */
    fn inv_top_left_matrix(&self) -> Matrix2<f64> {
        let m = self.matrix();
        Matrix2::new(
            m[(0, 0)], m[(0, 1)],
            m[(1, 0)], m[(1, 1)],
        ).try_inverse().unwrap()
    }

    /**
     * Top right 2x2 entries of M, defined in PipeSection::matrix(). Reused in several places.
     */
    fn top_right_matrix(&self) -> Matrix2<f64> {
        let m = self.matrix();
        Matrix2::new(
            m[(0, 2)], m[(0, 3)],
            m[(1, 2)], m[(1, 3)],
        )
    }

    /**
     * Inverse of PipeSection::matrix.
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
     * transformation T such that T(C) is this pipe section.
     */
    pub fn transformation_from_base_cylinder(&self) -> Affine3<f64> {
        Affine3::from_matrix_unchecked(self.inv_matrix())
    }

    /**
     * Define the base cylinder C = {(x, y, z): x^2 + y^2 = 1, z in R}. Return the affine
     * transformation T such that applying T to this pipe section produces C.
     */
    pub fn transformation_to_base_cylinder(&self) -> Affine3<f64> {
        Affine3::from_matrix_unchecked(self.matrix())
    }

    pub fn basic_forward_matrix(&self) -> Matrix3<f64> {
        let m = self.inv_matrix();
        Matrix3::new(
            m[(0, 0)], m[(0, 1)], m[(0, 2)],
            m[(1, 0)], m[(1, 1)], m[(1, 2)],
            m[(2, 0)], m[(2, 1)], m[(2, 2)],
        )
    }

    pub fn basic_backward_matrix(&self) -> Matrix3<f64> {
        let m = self.matrix();
        Matrix3::new(
            m[(0, 0)], m[(0, 1)], m[(0, 2)],
            m[(1, 0)], m[(1, 1)], m[(1, 2)],
            m[(2, 0)], m[(2, 1)], m[(2, 2)],
        )
    }

    /**
     * Compute the two displacement vectors that form the elliptic cross-section of the pipe.
     */
    fn ellipse_vertex_displacements(&self) -> (Vector3<f64>, Vector3<f64>) {
        let (_, direction) = self.axis_line();
        let z = Vector3::z();
        // R rotates the pipe's axis of symmetry to (0, 0, 1).
        let r = Rotation3::rotation_between(&direction, &z).unwrap_or(Rotation3::identity());
        let r_inv = Rotation3::rotation_between(&z, &direction).unwrap_or(Rotation3::identity());
        let project_xy: Matrix3<f64> = Matrix3::from_diagonal(&Vector3::new(1.0, 1.0, 0.0));
        let x = Vector3::x();
        let y = Vector3::y();
        let da = r_inv * project_xy * r * self.transformation_from_base_cylinder().transform_vector(&x);
        let db = r_inv * project_xy * r * self.transformation_from_base_cylinder().transform_vector(&y);
        if db.norm_squared() > da.norm_squared() {
            (db, da)
        } else {
            (da, db)
        }
    }

    fn ellipse_vertices(&self) -> (Vector3<f64>, Vector3<f64>) {
        let (start, _) = self.axis_line();
        let (da, db) = self.ellipse_vertex_displacements();
        (start + da, start + db)
    }

    fn surface_coords_to_cartesian(&self, u: f64, theta: f64) -> Vector3<f64> {
        let (start, direction) = self.axis_line();
        let (da, db) = self.ellipse_vertex_displacements();
        let a = da.norm();
        let b = db.norm();
        let theta_warped = warp_elliptic_angle(theta, a, b);
        let cx = theta_warped.cos();
        let cy = theta_warped.sin();
        start + u * direction + cx * da + cy * db
    }

    fn ellipse_circumference(&self) -> f64 {
        let (da, db) = self.ellipse_vertex_displacements();
        ellipse_circumference(da.norm(), db.norm())
    }

    pub fn as_mesh(&self) -> Mesh {
        if self.d != 0.0 {
            panic!("PipeSection::as_mesh does not yet work with d != 0");
        }

        if self.a == 0.0 && self.b == 0.0 {
            let tmp = 1.0 / squared(self.c) - squared(self.w);
            if tmp <= 0.0 {
                return Mesh::empty();
            }
            let z = tmp.sqrt();
            return Mesh::merge(vec![
                Mesh::plane(z),
                Mesh::plane(-z),
            ])
        }

        let spacing = 0.05;
        let radial_segments = self.ellipse_circumference() / spacing;
        let radial_segments = ((radial_segments / 4.0).ceil() * 4.0) as usize;

        let half_height = 2.0;
        let linear_segments = half_height * 2.0 / spacing;
        let linear_segments = ((linear_segments / 2.0).ceil() * 2.0) as usize;

        // Vertex indices: i * radial_segments + j
        let mut vertices: Vec<Vector3<f64>> = vec![];
        for i in 0..=linear_segments {
            let i_unipolar = (i as f64) / (linear_segments as f64);
            let i_bipolar = i_unipolar * 2.0 - 1.0;
            for j in 0..radial_segments {
                let theta = (j as f64) * std::f64::consts::TAU / (radial_segments as f64);
                let u = i_bipolar * half_height; 
                vertices.push(
                    self.surface_coords_to_cartesian(u, theta)
                );
            }
        }

        let mut faces: Vec<Face> = vec![];
        for i in 0..linear_segments {
            for j in 0..radial_segments {
                let v1 = i * radial_segments + j;
                let v2 = i * radial_segments + (j + 1) % radial_segments;
                let v3 = (i + 1) * radial_segments + j;
                let v4 = (i + 1) * radial_segments + (j + 1) % radial_segments;

                // v1 -- v2
                // | ,--' |
                // v3 -- v4
                faces.push(Face {
                    v1: v1,
                    v2: v2,
                    v3: v3,
                });
                faces.push(Face {
                    v1: v2,
                    v2: v4,
                    v3: v3,
                });
            }
        }
        Mesh { vertices, faces }
    }
}

#[cfg(test)]
mod test {
    use na::Matrix3;
    use super::*;

    #[test]
    fn test_pipe_axis() {
        let pipe = PipeSection { a: 1.2, b: -0.5, c: 0.4, d: 0.1, w: 0.1 };
        let (start, direction) = pipe.axis_line();
        assert_abs_diff_eq!(pipe.scalar_field(&start), -1.0);
        let point_2 = start + direction;
        assert_abs_diff_eq!(pipe.scalar_field(&point_2), -1.0);
        let point_3 = start + direction * 1.2345;
        assert_abs_diff_eq!(pipe.scalar_field(&point_3), -1.0);
    }

    #[test]
    fn test_inv_matrix() {
        let pipe = PipeSection { a: 1.2, b: -0.5, c: 0.4, d: 0.1, w: 0.1 };
        assert_abs_diff_eq!(
            pipe.matrix() * pipe.inv_matrix(),
            Matrix4::identity()
        )
    }

    #[test]
    fn test_ellipse_vertices() {
        let pipe = PipeSection { a: 1.2, b: -0.5, c: 0.4, d: 0.0, w: 0.1 };
        let (start, direction) = pipe.axis_line();
        let (da, db) = pipe.ellipse_vertex_displacements();
        let (pa, pb) = pipe.ellipse_vertices();
        assert_abs_diff_eq!(pipe.scalar_field(&pa), 0.0, epsilon = 1e-5);
        assert_abs_diff_eq!(pipe.scalar_field(&pb), 0.0, epsilon = 1e-5);
        assert_abs_diff_eq!(da.dot(&direction), 0.0, epsilon = 1e-5);
        assert_abs_diff_eq!(db.dot(&direction), 0.0, epsilon = 1e-5);
    }

    #[test]
    fn test_surface_to_cartesian() {
        let pipe = PipeSection { a: 1.2, b: -0.5, c: 0.4, d: 0.0, w: 0.1 };
        let point = pipe.surface_coords_to_cartesian(-1.34, 0.3);
        assert_abs_diff_eq!(pipe.scalar_field(&point), 0.0, epsilon = 1e-5);
    }
}