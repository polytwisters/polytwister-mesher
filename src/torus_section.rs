enum TorusSection {
    CylinderCylinder(Cylinder, Cylinder),
    CylinderZPlanes(Cylinder, f64),
    Empty,
}

impl TorusSection {
    fn new(pipe_section_1: PipeSection, pipe_section_2: PipeSection) -> Self {
        if pipe_section_1.is_plane() {
            if let Some(z) = pipe_section_1.plane_z() {
                if pipe_section_2.is_plane() {
                    // Pipe sections 1 and 2 are both planes.
                    TorusSection::Empty
                } else {
                    // Cylinder and z-planes intersection.
                    TorusSection::CylinderZPlanes(pipe_section_2.as_cylinder(), z)
                }
            } else {
                // Pipe section 1 is empty.
                TorusSection::Empty
            }
        } else {
            if pipe_section_2.is_plane() {
                if let Some(z) = pipe_section_2.plane_z() {
                    // Cylinder and z-planes intersection.
                    TorusSection::CylinderZPlanes(pipe_section_1.as_cylinder(), z)
                } else {
                    // Pipe section 2 is planar and empty.
                    TorusSection::Empty
                }
            } else {
                // Two cylinders.
                TorusSection::CylinderCylinder(
                    pipe_section_1.as_cylinder(),
                    pipe_section_2.as_cylinder()
                )
            }
        }
    }

    fn as_polylines(&self, resolution: usize) -> Vec<Polyline> {
        match self {
            TorusSection::Empty => vec![Polyline { points: vec![] }],
            TorusSection::CylinderCylinder(cylinder_1, cylinder_2) => {
                cylinder_1.intersect_cylinder(cylinder_2, resolution)
            },
            TorusSection::CylinderZPlanes(cylinder, z) => {
                cylinder.intersect_z_planes(*z, resolution)
            },
        }
    }
}