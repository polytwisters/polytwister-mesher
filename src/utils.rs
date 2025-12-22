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