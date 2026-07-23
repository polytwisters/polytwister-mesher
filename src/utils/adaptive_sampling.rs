#[derive(Clone, Copy, Debug)]
pub enum Topology1D {
    Linear,
    Circular
}

pub const TOLERANCE: f64 = 0.3;

pub struct AdaptiveSamplingConfig {
    pub max_distance: f64,
    pub distance_relative_tolerance: f64,
    pub t_range: (f64, f64),
    pub guess_num_points: usize,
    pub topology: Topology1D,
}

impl AdaptiveSamplingConfig {
    fn min_distance(&self) -> f64 {
        self.max_distance * (1.0 - self.distance_relative_tolerance)
    }
}

/**
 * Given a distance function d(t1, t2) and an interval on the real line, produce a list of t values
 * in that interval such that for all consecutive t1, t2 we have d_min <= d(t1, t2) <= d_max. The
 * list of t values may be linear (range interval is closed and the output contains both endpoints),
 * or circular (interval is half-open).
 */
pub fn adaptive_sample<F : Fn (f64, f64) -> f64>(
    distance_func: F,
    config: &AdaptiveSamplingConfig
) -> Vec<f64> {
    let (min, max) = config.t_range;
    let distance = config.max_distance;
    let min_distance = config.min_distance();
    let max_distance = config.max_distance;
    let topology = config.topology;
    let guess_num_points = config.guess_num_points;

    // Start with an evenly spaced number of points. In the circular case, they are evenly spaced
    // including min but excluding max. In the linear case, they are evenly spaced including both,
    // hence the conditional.
    let evenly_spaced: Vec<_> = (0..guess_num_points).map(|i| {
        let t = i as f64 / (match topology {
            Topology1D::Circular => guess_num_points,
            Topology1D::Linear => guess_num_points - 1,
        }) as f64;
        min + (max - min) * t
    }).collect();

    // Walk through each line segment connecting pairs of points and check distances.
    // Heuristic to reduce unnecessary reallocations: assume the final number of points needed is
    // about twice the initial num points.
    let mut result: Vec<f64> = Vec::with_capacity(evenly_spaced.len() * 2);
    for i1 in 0usize..evenly_spaced.len() {
        // Always add the point from the original.
        let x1 = evenly_spaced[i1];
        result.push(x1);

        let mut i2 = i1 + 1;
        // In the linear case, the final segment is ignored.
        if matches!(topology, Topology1D::Linear) && i2 >= evenly_spaced.len() {
            break;
        }
        i2 = i2.rem_euclid(evenly_spaced.len());

        let x2 = evenly_spaced[i2];
        // If the distance between successive points is larger than the minimum, subdivide it
        // into smaller segments.
        let d = distance_func(x1, x2);
        if d > max_distance {
            let subdivisions = (d / max_distance).ceil() as usize;
            // Start with 1 here, as we already added theta1.
            for i in 1..subdivisions {
                let t = i as f64 / subdivisions as f64;
                let x = x1 * (1.0 - t) + x2 * t;
                result.push(x);
            }
        }
    }
    result.shrink_to_fit();

    result
}

#[cfg(test)]
mod test {
    use nalgebra::Point2;
    use super::*;

    fn test_basic() {
        let func = |t: f64| {
            Point2::new(t.sin(), (t * 3.0).cos())
        };
        let distance_func = |t1, t2| {
            na::distance(&func(t1), &func(t2))
        };
        let distance = 0.01;
        let config = AdaptiveSamplingConfig {
            max_distance: distance,
            distance_relative_tolerance: 0.1,
            guess_num_points: 10,
            t_range: (0.3, 3.0),
            topology: Topology1D::Linear,
        };
        let t = adaptive_sample(distance_func, &config);
        for i in 0..t.len() - 1 {
            assert!(t[i] < t[i + 1]);
            assert!(distance_func(t[i], t[i + 1]) < config.max_distance);
        }
    }
}