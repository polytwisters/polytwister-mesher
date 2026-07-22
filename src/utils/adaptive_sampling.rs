/**
 * Given a distance function d(t1, t2) and an interval on the real line, produce a list of t values
 * in that interval such that d(t1, t2) is almost always less than a given maximum distance. The
 * list of t values may be linear (interval is closed and the output contains both endpoints), or
 * circular (interval is half-open).
 * 
 * This works by first evenly spacing initial_num_points in the interval, and for each pair of
 * consecutive t1, t2, if d(t1, t2) > max_distance then subdivide the interval [t1, t2] into
 * ceil(max_distance / d(t1, t2)) segments. This does not formally guarantee that the maximum
 * distance is adhered to but in practice it basically always works if the curve is reasonably
 * smooth and the number of initial points is sufficient.
 */
pub fn adaptive_sample<F : Fn (f64, f64) -> f64>(
    distance_func: F,
    max_distance: f64,
    initial_num_points: usize,
    min: f64,
    max: f64,
    circular: bool
) -> Vec<f64> {
    // Start with an evenly spaced number of points. In the circular case, they are evenly spaced
    // including min but excluding max. In the linear case, they are evenly spaced including both,
    // hence the conditional.
    let evenly_spaced: Vec<_> = (0..initial_num_points).map(|i| {
        let t = i as f64 / (if circular { initial_num_points } else { initial_num_points - 1 }) as f64;
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
        if !circular && i2 >= evenly_spaced.len() {
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
    use crate::utils::adaptive_sample;

    fn test_basic() {
        let func = |x: f64| { x.sin() };
        let distance_func = |t1, t2| {
            let p1 = Point2::new(t1, func(t2));
            let p2 = Point2::new(t2, func(t2));
            na::distance(&p1, &p2)
        };
        let max_distance = 0.01;
        let t = adaptive_sample(
            distance_func,
            max_distance,
            10,
            0.3,
            3.0,
            false
        );
        for i in 0..t.len() - 1 {
            assert!(t[i] < t[i + 1]);
            assert!(distance_func(t[i], t[i + 1]) < max_distance);
        }
    }
}