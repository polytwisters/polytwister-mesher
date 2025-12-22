extern crate nalgebra as na;
use na::{Vector3, Matrix3, Rotation3};
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
    pub w: f64,
}


impl PipeSection {
    /**
     * Evaluate the scalar field associated with this pipe cross section:
     * 
     *     F(x, y, z) = (ax + by + cz)^2 + (bx - ay - cw)^2 - 1
     * 
     * The pipe cross section is given by the isosurface F(x, y, z) = 0. For other points, if
     * F(x, y, z) < 0 then the point is inside the pipe, F(x, y, z) > 0 is outside the pipe, and
     * F(x, y, z) = -1 is on the pipe's symmetry axis (plane of symmetry if a = b = 0).
     */
    pub fn scalar_field(&self, point: &Vector3<f64>) -> f64 {
        squared(self.a * point.x + self.b * point.y + self.c * point.z)
        + squared(self.b * point.x - self.a * point.y - self.c * self.w)
        - 1.0
    }

    fn abcw(&self) -> (f64, f64, f64, f64) {
        (self.a, self.b, self.c, self.w)
    }

    /**
     * Return an explicit parametrization of the line which is this pipe section's symmetry axis.
     * The parametrization is v(t) = v_0 + d * t and returned as (v_0, d) so that v_0 is the
     * starting point and d is the direction vector. d is always a unit vector.
     */
    fn axis_line(&self) -> (Vector3<f64>, Vector3<f64>) {
        let (a, b, c, w) = self.abcw();
        let tmp = 1.0 / (squared(a) + squared(b));
        let direction_vector = Vector3::new( -a * c * tmp, -b * c * tmp, 1.0).normalize();
        let start = Vector3::new(b * c * w * tmp, -a * c * w * tmp, 0.0);
        return (start, direction_vector);
    }

    /**
     * Return a 3x3 matrix that turns this pipe into the "base cylinder" x^2 + y^2 = 1, z in R.
     * It is assumed that a = b = 0 does not hold. Translation due to w is ignored.
     * 
     * This is called the "basic" backward matrix because it comes directly from the implicit
     * equations. Although it transforms the z-axis to the cylinder's axis of symmetry, it does not
     * preserve distance along that line.
     */
    fn basic_backward_matrix(&self) -> Matrix3<f64> {
        let (a, b, c, _) = self.abcw();
        Matrix3::new(
            a, b, c,
            b, -a, 0.0,
            0.0, 0.0, 1.0,
        )
    }

    /**
     * Return a 3x3 matrix that turns the "base cylinder" x^2 + y^2 = 1, z in R into this pipe.
     * It is assumed that a = b = 0 does not hold. Translation due to w is ignored.
     * 
     * This is the inverse of the basic_backward_matrix.
     */
    pub fn basic_forward_matrix(&self) -> Matrix3<f64> {
        let (a, b, c, _) = self.abcw();
        let tmp = 1.0 / (squared(a) + squared(b));
        Matrix3::new(
            a * tmp, b * tmp, -a * c * tmp,
            b * tmp, -a * tmp, -b * c * tmp,
            0.0, 0.0, 1.0,
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
        let da = r_inv * project_xy * r * self.basic_forward_matrix() * Vector3::x();
        let db = r_inv * project_xy * r * self.basic_forward_matrix() * Vector3::y();
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
        let pipe = PipeSection { a: 1.2, b: -0.5, c: 0.4, w: 0.1 };
        let (start, direction) = pipe.axis_line();
        assert_abs_diff_eq!(pipe.scalar_field(&start), -1.0);
        let point_2 = start + direction;
        assert_abs_diff_eq!(pipe.scalar_field(&point_2), -1.0);
        let point_3 = start + direction * 1.2345;
        assert_abs_diff_eq!(pipe.scalar_field(&point_3), -1.0);
    }

    #[test]
    fn test_forward_inverse_matrices() {
        let pipe = PipeSection { a: 1.2, b: -0.5, c: 0.4, w: 0.1 };
        assert_abs_diff_eq!(
            pipe.basic_forward_matrix() * pipe.basic_backward_matrix(),
            Matrix3::identity()
        )
    }

    #[test]
    fn test_ellipse_vertices() {
        let pipe = PipeSection { a: 1.2, b: -0.5, c: 0.4, w: 0.1 };
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
        let pipe = PipeSection { a: 1.2, b: -0.5, c: 0.4, w: 0.1 };
        let point = pipe.surface_coords_to_cartesian(-1.34, 0.3);
        assert_abs_diff_eq!(pipe.scalar_field(&point), 0.0, epsilon = 1e-5);
    }
}