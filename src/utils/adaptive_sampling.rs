use crate::utils::{bisection_search, lerp, lerp_inverse, linspace};

pub enum Topology1D {
    Circular,
    Linear
}

fn binary_search(array: &Vec<f32>, target: f32) -> usize {
    let len = array.len();
    assert!(len >= 2);
    let mut min = 0;
    let mut max = array.len() - 1;
    loop {
        if target < array[min] {
            return min;
        }
        if target >= array[max] {
            return max.min(len - 2);
        }
        if !(min + 1 < max) {
            break;
        }
        let midpoint_index = min + (max - min) / 2;
        if midpoint_index == min || midpoint_index == max {
            break;
        }
        let midpoint_value = array[midpoint_index];
        if midpoint_value <= target {
            min = midpoint_index;
        } else {
            max = midpoint_index;
        }
    }
    min
}

fn linear_interpolate(t: &Vec<f32>, f: &Vec<f32>, t_in: f32) -> f32 {
    let i1 = binary_search(t, t_in);
    let i2 = i1 + 1;
    let (t1, t2) = (t[i1], t[i2]);
    let (f1, f2) = (f[i1], f[i2]);
    let frac = lerp_inverse(t1, t2, t_in);
    lerp(f1, f2, frac)
}

pub struct AdaptiveSamplingConfig {
    pub target_distance: f32,
    pub t_range: (f32, f32),
    pub guess_num_points: usize,
}

/**
 * Given a distance function d(t1, t2), produce a list of t values such that:
 * 
 * * The first value is t_range.0.
 * * All values are in increasing order.
 * * If the topology is linear, then the final value is t_range.1.
 * * The values are roughly evenly spaced and their distance is roughly target_distance.
 * * If the topology is circular, all values are less than t_range.1, and the distance between the
 * last and first points is accounted for.
 * 
 * If the topology is circular then it is assumed that d(t_range.0, t_range.1) is 0.
 */
pub fn adaptive_sample<F : Fn (f32, f32) -> f32>(
    distance_func: F,
    topology: Topology1D,
    config: &AdaptiveSamplingConfig,
) -> Vec<f32> {
    let (min, max) = config.t_range;
    let max_distance = config.target_distance;
    let guess_num_points = config.guess_num_points;

    // Start with an evenly spaced number of points. Even in the circular case we want to include
    // the max point.
    let evenly_spaced_ts = linspace(min, max, guess_num_points);

    // Compute distances between consecutive points and take the cumulative sum of them to produce
    // the full arc length.
    let (cumulative_arc_lengths, total_arc_length) = {
        let mut result = Vec::with_capacity(evenly_spaced_ts.len());
        let mut cumulative_arc_length = 0.0;
        for i1 in 0..evenly_spaced_ts.len() - 1 {
            let i2 = i1 + 1;
            let arc_length = distance_func(evenly_spaced_ts[i1], evenly_spaced_ts[i2]);
            result.push(cumulative_arc_length);
            cumulative_arc_length += arc_length;
        }
        (result, cumulative_arc_length)
    };

    // There needs to be at least one segment. Three is a nicer minimum because for a closed 3D
    // polyline, you want at least three points.
    let num_segments = (
        (total_arc_length / max_distance).ceil() as usize
    ).max(3);
    let num_points = match topology {
        Topology1D::Circular => num_segments,
        Topology1D::Linear => num_segments + 1,
    };

    let result: Vec<_> = (0..num_points).map(|j| {
        let arc_length = j as f32 / num_segments as f32 * total_arc_length;
        linear_interpolate(&cumulative_arc_lengths, &evenly_spaced_ts, arc_length)
    }).collect();

    result
}


#[cfg(test)]
mod test {
    use core::f32;
use std::cmp::max;

use nalgebra::Point2;
    use super::*;

    /**
     * Using an arbitrary 2D plane curve as an example, verify that the points are in range, are
     * close to each other, and points with spacing close to the target distance.
     */
    #[test]
    fn test_linear() {
        let func = |t: f32| {
            Point2::new(t.sin(), (t * 3.0).cos())
        };
        let distance_func = |t1, t2| {
            na::distance(&func(t1), &func(t2))
        };
        let target_distance = 0.1;
        let t_range = (0.3, 3.0);
        let config = AdaptiveSamplingConfig {
            target_distance,
            guess_num_points: 100,
            t_range,
        };
        let t = adaptive_sample(distance_func, Topology1D::Linear, &config);
        let mut min_distance = f32::INFINITY;
        let mut max_distance = f32::NEG_INFINITY;
        assert!(t.len() > 2);
        assert!(t.iter().all(|&t| t_range.0 <= t && t < t_range.1));
        for i in 0..t.len() - 1 {
            assert!(t[i] < t[i + 1]);
            let distance = distance_func(t[i], t[i + 1]);
            min_distance = distance.min(min_distance);
            max_distance = distance.max(max_distance);
        }
        let min_distance_ratio = min_distance / target_distance;
        let max_distance_ratio = max_distance / target_distance;
        let relative_window_size = 0.5;
        assert!(min_distance_ratio > 1.0 - relative_window_size);
        assert!(max_distance_ratio < 1.0 + relative_window_size);
    }

    /// If the overall arc length is very short, make sure at least two points are returned.
    #[test]
    fn test_very_short() {
        let distance_func = |t1: f32, t2: f32| {
            (t2 - t1).abs()
        };
        let target_distance = 0.1;
        let t_range = (0.0, 0.01);
        let config = AdaptiveSamplingConfig {
            target_distance,
            guess_num_points: 100,
            t_range,
        };
        let t = adaptive_sample(distance_func, Topology1D::Linear, &config);
        assert!(t.len() > 2);
    }
}