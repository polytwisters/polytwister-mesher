#[derive(Clone, Debug)]
pub struct Grid {
    pub u_axis: GridAxis,
    pub v_axis: GridAxis,
}

#[derive(Clone, Copy, Debug)]
pub enum GridAxisTopology {
    Linear,
    Circular
}

/**
 * A single axis of a one-dimensional grid, with support for non-uniform grid sampling and circular
 * or linear topologies.
 * 
 * The grid axis's function is to convert from "grid coordinates" to "values." For example, a size-3
 * uniform grid might have grid coordinates [0, 1, 2] mapped to values [0.0, 0.5, 1.0].
 */
#[derive(Clone, Debug)]
pub struct GridAxis {
    topology: GridAxisTopology,
    size: usize,
    values: Vec<f64>,
}

impl GridAxis {
    pub fn uniform_linear(size: usize, min: f64, max: f64) -> Self {
        Self {
            topology: GridAxisTopology::Linear,
            size,
            values: (0..size).map(|i| {
                let t = (i as f64) / (size - 1) as f64;
                min + t * (max - min)
            }).collect()
        }
    }

    pub fn uniform_circular(size: usize, max: f64) -> Self {
        Self {
            topology: GridAxisTopology::Circular,
            size,
            values: (0..size).map(|i| (i as f64) / size as f64 * max).collect()
        }
    }

    /// Convert from grid coordinate to values. If it is on an integer grid point, return the value
    /// at that grid point. If it is between two grid points, linearly interpolate the values.
    pub fn at(&self, index: f64) -> f64 {
        match self.topology {
            GridAxisTopology::Circular => {
                let index1 = self.wrap(index as usize);
                let index2 = self.wrap(index as usize + 1);
                let t = index.fract();
                self.values[index1] * (1.0 - t) + self.values[index2] * t
            },
            GridAxisTopology::Linear => {
                let mut index1 = (index as usize).min(self.size - 2);
                let index2 = index1 + 1;
                let t = index - index1 as f64;
                self.values[index1] * (1.0 - t) + self.values[index2] * t
            },
        }
    }

    pub fn num_segments(&self) -> usize {
        match self.topology {
            GridAxisTopology::Circular => self.size,
            GridAxisTopology::Linear => self.size - 1,
        }
    }

    /// For a linear GridAxis, do nothing. For a circular grid axis, take the index modulo the size
    /// of the grid.
    pub fn wrap(&self, index: usize) -> usize {
        match self.topology {
            GridAxisTopology::Circular => index.rem_euclid(self.size),
            GridAxisTopology::Linear => index,
        }
    }
}