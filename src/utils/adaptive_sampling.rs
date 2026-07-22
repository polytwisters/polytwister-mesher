pub fn adaptive_sample<F : Fn (f64, f64) -> f64>(
    distance_func: F,
    max_distance: f64,
    initial_num_points: usize,
    min: f64,
    max: f64,
    circular: bool
) -> Vec<f64> {
    let initial: Vec<_> = (0..initial_num_points).map(|i| {
        let t = i as f64 / (if circular { initial_num_points } else { initial_num_points - 1 }) as f64;
        min + (max - min) * t
    }).collect();

    // Walk through each line segment connecting pairs of points and check distances.

    // Rough heuristic, assume the number of new points needed is about twice.
    let mut result: Vec<f64> = Vec::with_capacity(initial.len() * 2);
    for i1 in 0usize..initial.len() {
        let x1 = initial[i1];
        // Always keep the first point in the segment.
        result.push(x1);

        let mut i2 = i1 + 1;
        if circular && i2 >= initial.len() {
            break;
        }
        i2 = i2.rem_euclid(initial.len());

        let x2 = initial[i2];
        let d = distance_func(x1, x2);
        // If the distance between successive points is larger than the minimum, subdivide it
        // into smaller segments.
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