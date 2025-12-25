extern crate nalgebra as na;
use na::{Affine3, Matrix2, Matrix3, Matrix4, Rotation3, Vector2, Vector3, Point3};
use crate::cylinder::{Cylinder, CylinderMeshOptions};
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
    pub fn new(a: f64, b: f64, c: f64, d: f64, w: f64) -> Self {
        PipeSection { a, b, c, d, w }
    }

    /**
     * Evaluate the scalar field associated with this pipe cross section:
     * 
     *     F(x, y, z) = (ax + by + cz + dw)^2 + (bx - ay + dz - cw)^2 - 1
     * 
     * The pipe cross section is given by the isosurface F(x, y, z) = 0. For other points, if
     * F(x, y, z) < 0 then the point is inside the pipe, F(x, y, z) > 0 is outside the pipe, and
     * F(x, y, z) = -1 is on the pipe's symmetry axis (plane of symmetry if a = b = 0).
     */
    pub fn scalar_field(&self, point: &Point3<f64>) -> f64 {
        squared(self.a * point.x + self.b * point.y + self.c * point.z + self.d * self.w)
        + squared(self.b * point.x - self.a * point.y + self.d * point.z - self.c * self.w)
        - 1.0
    }

    pub fn contains(&self, point: &Point3<f64>) -> bool {
        self.scalar_field(point) < 0.0
    }

    pub fn as_cylinder(&self) -> Cylinder {
        Cylinder {
            m11: self.a,
            m12: self.b,
            m13: self.c,
            m14: self.d * self.w,
            m21: self.b,
            m22: -self.a,
            m23: self.d,
            m24: -self.c * self.w,
        }
    }

    pub fn as_mesh(&self, options: &CylinderMeshOptions) -> Mesh {
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
                Mesh::plane(z, options.half_length, options.linear_segments),
                Mesh::plane(-z, options.half_length, options.linear_segments),
            ])
        }

        self.as_cylinder().as_mesh(options)
    }
}