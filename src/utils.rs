use std::f64;
use nalgebra::{Point2};

pub fn squared(x: f64) -> f64 {
    x * x
}

pub fn angle(point: &Point2<f64>) -> f64 {
    f64::atan2(point.y, point.x).rem_euclid(f64::consts::TAU)
}

pub fn sort2(x: (f64, f64)) -> (f64, f64) {
    let (x1, x2) = x;
    if x1 < x2 { (x1, x2) } else { (x2, x1) }
}

pub fn sort4(x: (f64, f64, f64, f64)) -> (f64, f64, f64, f64) {
    let mut tmp = [x.0, x.1, x.2, x.3];
    tmp.sort_by(f64::total_cmp);
    (tmp[0], tmp[1], tmp[2], tmp[3])
}

/**
 * Return a linearly spaced series of n values in the closed interval [start, end]. Both endpoints
 * are inclusive.
 */
pub fn linspace(start: f64, end: f64, n: usize) -> Vec<f64> {
    (0..n).map(|i| start + (i as f64) / ((n - 1) as f64) * (end - start)).collect::<_>()
}

/**
 * Convert the vector of pairs: [(1, X), (2, 10), ..., (4, 6), (X, 5)] to the vector [1, 2, 3, ...].
 */
pub fn unzip_circle<T>(pairs: Vec<(T, T)>) -> Vec<T> {
    let (mut tmp1, mut tmp2): (Vec<T>, Vec<T>) = pairs.into_iter().unzip();
    tmp1.pop();
    tmp2.reverse();
    tmp2.pop();
    tmp1.append(&mut tmp2);
    tmp1
}