extern crate nalgebra as na;
use core::f64;

use na::{Affine3, Matrix2, Matrix3, Matrix4, Point2, Point3, Rotation3, Vector2, Vector3};
use crate::pipe_section::PipeSection;
use crate::utils::{Ellipse, squared};
use crate::utils::adaptive_sampling::{self, AdaptiveSamplingConfig, adaptive_sample};
use crate::mesh::{Face, Mesh, Vertex};
use crate::config::{CylinderMeshConfig};

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
            m11: 0.3,
            m12: -0.2,
            m13: 0.4,
            m14: -0.3,
            m21: 0.2,
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
     * Gradient of the scalar field.
     */
    pub fn scalar_field_gradient(&self, point: &Point3<f64>) -> Vector3<f64> {
        let tmp1 = self.m11 * point.x + self.m12 * point.y + self.m13 * point.z + self.m14;
        let tmp2 = self.m21 * point.x + self.m22 * point.y + self.m23 * point.z - self.m24;
        Vector3::new(
            2.0 * tmp1 * self.m11 + 2.0 * tmp2 * self.m21,
            2.0 * tmp1 * self.m12 + 2.0 * tmp2 * self.m22,
            2.0 * tmp1 * self.m13 + 2.0 * tmp2 * self.m23,
        )
    }

    /// Return true if the point is on the boundary or interior of this cylinder.
    pub fn interior_contains(&self, point: &Point3<f64>) -> bool {
        self.scalar_field(point) <= 0.0
    }

    pub fn intersect_axis_line_z_plane(&self, z: f64) -> Point3<f64> {
        // Solve M(x, y, z, 1) = 0 with z fixed. Ignore bottom two rows and use block matrix
        // inversion:
        //     (x, y) = -inv_top_left * top_right (z, 1)
        let point_2d = -self.inv_top_left_matrix() * (self.top_right_matrix() * Vector2::new(z, 1.0));
        Point3::new(point_2d.x, point_2d.y, z)
    }

    /**
     * Return an explicit parametrization of the line which is this cylinder section's symmetry axis.
     * The parametrization is v(t) = v_0 + d * t and returned as (v_0, d) so that v_0 is the
     * starting point and d is the direction vector. d is always a unit vector.
     * 
     * This method assumes that the cylinder intersects the plane z = 0. This holds for all pipe
     * cross sections but not necessarily all affine transformations.
     */
    pub fn axis_line(&self) -> (Point3<f64>, Vector3<f64>) {
        let point1 = self.intersect_axis_line_z_plane(0.0);
        let point2 = self.intersect_axis_line_z_plane(1.0);
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
     * Top left 2x2 entries of M.
     */
    pub fn top_left_matrix(&self) -> Matrix2<f64> {
        let m = self.matrix();
        Matrix2::new(
            m[(0, 0)], m[(0, 1)],
            m[(1, 0)], m[(1, 1)],
        )
    }

    /**
     * The inverse of the top left 2x2 entries of M, defined in Cylinder::matrix(). Reused in
     * several places.
     */
    pub fn inv_top_left_matrix(&self) -> Matrix2<f64> {
        self.top_left_matrix().try_inverse().unwrap()
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
     * The inverse of the top left 3x3 entries of M, defined in Cylinder::matrix().
     */
    pub fn top_left_matrix_3(&self) -> Matrix3<f64> {
        let m = self.matrix();
        // Sorry about repetitive code -- I can NOT figure out how to take a static slice of a
        // Matrix4.
        Matrix3::new(
            m[(0, 0)], m[(0, 1)], m[(0, 2)],
            m[(1, 0)], m[(1, 1)], m[(1, 2)],
            m[(2, 0)], m[(2, 1)], m[(2, 2)],
        )
    }

    /**
     * The intersection of the cylinder with a plane orthogonal to the cylinder's axis of symmetry
     * is an ellipse. The vertices of an ellipse are the points furthest and closest to its center,
     * if they exist. Return the displacement vectors, in the 3D world coordinate space, from the
     * ellipse's center to its vertices. They are returned in the order
     * (major vertex displacement, minor vertex displacement). These two vectors are guaranteed
     * orthogonal to each other, and their cross product is in the same direction as the axis line.
     * 
     * If the ellipse is a circle, there are no defined vertices, so this method returns any two
     * vectors whose lengths are equal to the radius of the circle, and which are orthogonal to each
     * other and to the axis.
     */
    pub fn ellipse_vertex_displacements(&self) -> (Vector3<f64>, Vector3<f64>) {
        let (_, direction) = self.axis_line();
        // 3D rotation turning the z-axis into the axial line.
        let rotation = Rotation3::rotation_between(&Vector3::z(), &direction).unwrap_or(Rotation3::identity());
        // Rotation to the base cylinder. Translation is ignored.
        let t: Matrix3<f64> = self.top_left_matrix_3() * rotation;
        // The top left 2x2 of the matrix gives the coefficient of an ellipse.
        let ellipse = Ellipse {
            matrix: Matrix2::new(
                t[(0, 0)], t[(0, 1)],
                t[(1, 0)], t[(1, 1)],
            )
        };
        let (major_2d, minor_2d) = ellipse.vertices();
        let major = Vector3::new(major_2d.x, major_2d.y, 0.0);
        let minor = Vector3::new(minor_2d.x, minor_2d.y, 0.0);
        let (major, minor) = (rotation * major, rotation * minor);
        if major.cross(&minor).dot(&direction) > 0.0 {
            (major, -minor)
        } else {
            (major, minor)
        }
    }

    pub fn ellipse_radii(&self) -> (f64, f64) {
        let (major, minor) = self.ellipse_vertex_displacements();
        (major.norm(), minor.norm())
    }

    /**
     * Find the vertices of a cross-sectional ellipse.
     */
    fn ellipse_vertices(&self) -> (Point3<f64>, Point3<f64>) {
        let (start, _) = self.axis_line();
        let (da, db) = self.ellipse_vertex_displacements();
        (start + da, start + db)
    }

    /// Convert 2D surface coordinates (u, theta) on the cylinder to points in 3D space. Changing u
    /// by a factor of du moves the point parallel to the cylinder's axis by a distance of du, and
    /// the range of u is the entire real line. Changing theta moves the point in an elliptical loop
    /// orthogonal to the cylinder's axis, ranging from 0 to 2pi.
    pub fn surface_coords_to_cartesian(&self, u: f64, theta: f64) -> Point3<f64> {
        let (start, direction) = self.axis_line();
        let (da, db) = self.ellipse_vertex_displacements();
        let cx = theta.cos();
        let cy = theta.sin();
        start + u * direction + cx * da + cy * db
    }

    /// Convert Cartesian coordinates to cylindrical coordinates (u, theta, r). If r = 1 then this
    /// is the inverse of Cylinder::surface_coords_to_cartesian.
    pub fn cartesian_to_cylindrical(&self, point: &Point3<f64>) -> (f64, f64, f64) {
        let (start, direction_normalized) = self.axis_line();
        let (da, db) = self.ellipse_vertex_displacements();
        let tmp = point - start;
        let u = tmp.dot(&direction_normalized);
        let cos_theta = tmp.dot(&da) / da.norm_squared();
        let sin_theta = tmp.dot(&db) / db.norm_squared();
        let r = cos_theta.hypot(sin_theta);
        let theta = sin_theta.atan2(cos_theta).rem_euclid(f64::consts::TAU);
        (u, theta, r)
    }


    /// Discretize this cylinder as a mesh.
    pub fn as_mesh(&self, options: &CylinderMeshConfig) -> Mesh {
        self.as_mesh_partial(|_| true, &options)
    }

    /// Discretize this cylinder as a mesh, but only include the points for
    /// which the predicate returns true.
    pub fn as_mesh_partial<F: Fn(&Point3<f64>) -> bool>(
        &self,
        predicate: F,
        config: &CylinderMeshConfig
    ) -> Mesh {
        let resolution = config.resolution;
        let half_height = config.half_length;
        let linear_segments = (half_height * 2.0 / resolution).ceil() as usize;

        let thetas = self.sample_theta(resolution);
        let radial_segments = thetas.len();

        // Vertex indices: i * radial_segments + j
        let mut vertices: Vec<Option<Vertex>> = vec![];
        for i in 0..=linear_segments {
            let i_unipolar = (i as f64) / (linear_segments as f64);
            let i_bipolar = i_unipolar * 2.0 - 1.0;
            for j in 0..radial_segments {
                let theta = thetas[j];
                let u = i_bipolar * half_height; 
                let point = self.surface_coords_to_cartesian(u, theta);
                vertices.push(if predicate(&point) {
                    let normal = self.scalar_field_gradient(&point).normalize();
                    Some(Vertex::new(point, normal))
                } else { None });
            }
        }

        let mut faces: Vec<Face> = vec![];
        for i in 0..linear_segments {
            for j in 0..radial_segments {
                let v1 = i * radial_segments + j;
                let v2 = i * radial_segments + (j + 1) % radial_segments;
                let v3 = (i + 1) * radial_segments + j;
                let v4 = (i + 1) * radial_segments + (j + 1) % radial_segments;

                let v2_exists = matches!(vertices[v2], Some(_));
                let v3_exists = matches!(vertices[v3], Some(_));

                // If all 4 vertices are present, then we triangulate the square
                // with 2 triangles, each clockwise in this diagram:
                //
                // v1 -- v2
                // | ,--' |
                // v3 -- v4
                //
                // If v2 or v3 are not present, we use the alternate
                // triangulation:
                //
                // v1----v2
                // | `--. |
                // v3 -- v4
                //
                // Mesh::from_partial will delete all triangles with missing
                // vertices.

                if v2_exists && v3_exists {
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
                } else {
                    faces.push(Face {
                        v1: v1,
                        v2: v2,
                        v3: v4,
                    });
                    faces.push(Face {
                        v1: v1,
                        v2: v4,
                        v3: v3,
                    });
                }
            }
        }
        Mesh::from_partial(&vertices, &faces)
    }

    /**
     * Perform nonuniform sampling of theta values, creating points spaced around an ellipse so that
     * no two consecutive points are closer together than the resolution.
     */
    pub fn sample_theta(&self, resolution: f64) -> Vec<f64> {
        // I initially played with methods using inverse incomplete elliptic integrals to space
        // points on an ellipse, but it was too annoying to fiddle with the dependencies needed to
        // do that. This is a simple brute-force numerical solution.

        // Evenly space points on a circle whose radius is the average fo the major and minor radii.
        // Just a heuristic.
        let (major_radius, minor_radius) = self.ellipse_radii();
        let average_radius = (major_radius + minor_radius) / 2.0;
        let guess_num_points = (
            (average_radius * f64::consts::TAU / resolution as f64) as usize
        ).max(4);

        let result = adaptive_sample(
            |theta1, theta2| {
                let p1 = Point2::new(major_radius * theta1.cos(), minor_radius * theta1.sin());
                let p2 = Point2::new(major_radius * theta2.cos(), minor_radius * theta2.sin());
                na::distance(&p1, &p2)
            },
            &AdaptiveSamplingConfig {
                max_distance: resolution,
                distance_relative_tolerance: adaptive_sampling::TOLERANCE,
                t_range: (0.0, f64::consts::TAU),
                guess_num_points,
                topology: adaptive_sampling::Topology1D::Circular,
            },
        );
        result
    }
}

#[cfg(test)]
mod test {
    use crate::cylinder;

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

    #[test]
    fn test_basis() {
        let cylinder = example_cylinder();
        let (_, d) = cylinder.axis_line();
        let (da, db) = cylinder.ellipse_vertex_displacements();
        assert_abs_diff_eq!(d.dot(&da), 0.0);
        assert_abs_diff_eq!(d.dot(&db), 0.0);
        assert_abs_diff_eq!(da.dot(&db), 0.0);
        assert!(da.cross(&db).dot(&d) < 0.0);
    }

    #[test]
    fn test_surface_coords() {
        let cylinder = example_cylinder();
        let u = 0.75;
        let theta = 3.345;
        let point = cylinder.surface_coords_to_cartesian(u, theta);
        let (u_out, theta_out, r_out) = cylinder.cartesian_to_cylindrical(&point);
        assert_abs_diff_eq!(u, u_out, epsilon = 1e-10);
        assert_abs_diff_eq!(r_out, 1.0, epsilon = 1e-10);
        assert_abs_diff_eq!(theta, theta_out, epsilon = 1e-10);
    }

    #[test]
    fn test_sample_theta() {
        let cylinder = example_cylinder();
        let min_distance = 0.05;
        let thetas = cylinder.sample_theta(min_distance);
        let u = 0.0;  // Doesn't matter.
        for i1 in 0..thetas.len() {
            let i2 = (i1 + 1).rem_euclid(thetas.len());
            let p1 = cylinder.surface_coords_to_cartesian(u, thetas[i1]);
            let p2 = cylinder.surface_coords_to_cartesian(u, thetas[i2]);
            assert!(na::distance(&p1, &p2) <= min_distance);
        }
    }
}