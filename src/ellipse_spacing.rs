use core::{f64};
use crate::utils::squared;

fn k_to_blend(k: f64) -> f64 {
    let m = k * k;
    let exponent = 0.7224221;
    let fade = 2.39077584;
    (1.0 - (1.0 - m).powf(exponent)) * fade + m * (1.0 - fade)
}

fn ellipse_unwarp_core(x: f64, k: f64) -> f64 {
    let b = k_to_blend(k);
    (1.0 - b) * x + b * (x * f64::consts::FRAC_PI_2).sin()
}

fn ellipse_unwarp_core_derivative(x: f64, k: f64) -> f64 {
    let b = k_to_blend(k);
    (1.0 - b) + b * (x * f64::consts::FRAC_PI_2).cos() * f64::consts::FRAC_PI_2
}

/**
 * Map the interval [0, 1] -> [0, 1] with a warping function. Given an ellipse parametrized as
 * (x, y) = (a sin theta, b cos theta) with a > b, with eccentricity k = sqrt(1 - b^2/a^2),
 * let theta = warp(i / (pi / 2)) * (pi / 2) with i ranging from 0 to 1. If i values are evenly
 * spaced, the points on the ellipse are evenly spaced.
 */
fn ellipse_warp_core(y: f64, k: f64) -> f64 {
    let mut x = y;
    for _ in 0..4 {
        let error = ellipse_unwarp_core(x, k) - y;
        if error.abs() < 1e-5 {
            return x;
        }
        let derivative_error = ellipse_unwarp_core_derivative(x, k);
        x -= error / derivative_error;
    }
    x
}

fn ellipse_warp_core_flip(q: f64, k: f64) -> f64 {
    1.0 - ellipse_warp_core(1.0 - q, k)
}

// https://www.e-magnetica.pl/doku.php/approximation_of_complete_elliptic_integrals
fn elliptic_e_complete(k: f64) -> f64 {
    f64::consts::FRAC_PI_2 - 0.567 * k.powf(2.4 + (k + 0.1).powf(5.8))
}

pub fn ellipse_circumference(width: f64, height: f64) -> f64 {
    let (a, b) = if width >= height { (width, height) } else { (height, width) };
    let k = (1.0 - squared(b / a)).sqrt();
    4.0 * a * elliptic_e_complete(k)
}

pub fn warp_elliptic_angle(phi: f64, a: f64, b: f64) -> f64 {
    // q = number of quarter turns
    let q = phi / f64::consts::FRAC_PI_2;
    let qw = if a >= b {
        let k = (1.0 - squared(b / a)).sqrt();
        match q as u8 {
            0 => ellipse_warp_core_flip(q, k),
            1 => 1.0 + ellipse_warp_core(q - 1.0, k),
            2 => 2.0 + ellipse_warp_core_flip(q - 2.0, k),
            _ => 3.0 + ellipse_warp_core(q - 3.0, k),
        }
    } else {
        let k = (1.0 - squared(a / b)).sqrt();
        match q as u8 {
            0 => ellipse_warp_core(q, k),
            1 => 1.0 + ellipse_warp_core_flip(q - 1.0, k),
            2 => 2.0 + ellipse_warp_core(q - 2.0, k),
            _ => 3.0 + ellipse_warp_core_flip(q - 3.0, k),
        }
    };
    qw * f64::consts::FRAC_PI_2
}

fn evenly_spaced_ellipse_points(a: f64, b: f64, n: usize) -> Vec<(f64, f64)> {
    (0..n).into_iter().map(|i| {
        let phi = (i as f64) / (n as f64) * f64::consts::TAU;
        let phi2 = warp_elliptic_angle(phi, a, b);
        (a * phi2.cos(), b * phi2.sin())
    }).collect::<Vec<(f64, f64)>>()
}

fn naively_spaced_ellipse_points(a: f64, b: f64, n: usize) -> Vec<(f64, f64)> {
    (0..n).into_iter().map(|i| {
        let phi = (i as f64) / (n as f64) * f64::consts::TAU;
        (a * phi.cos(), b * phi.sin())
    }).collect::<Vec<(f64, f64)>>()
}

#[cfg(test)]
mod test {
    use super::*;
    use std::iter::zip;

    #[test]
    fn test_ellip_inverse() {
        let k = 0.95;
        let x = 0.45;
        let e = ellipse_warp_core(x, k);
        let x2 = ellipse_unwarp_core(e, k);
        assert_abs_diff_eq!(x, x2, epsilon = 1e-5);
    }

    fn consecutive_distances(points: &Vec<(f64, f64)>) -> Vec<f64> {
        zip(&points[1..], &points[..points.len() - 1]).map(|((x1, y1), (x2, y2))| {
            f64::hypot(x1 - x2, y1 - y2)
        }).collect::<Vec<f64>>()
    }

    fn distance_range(points: &Vec<(f64, f64)>) -> f64 {
        let distances = consecutive_distances(&points);
        distances.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b))
        - distances.iter().fold(f64::INFINITY, |a, &b| a.min(b))
    }

    #[test]
    fn test_evenly_spaced_points() {
        let a = 3.0;
        let b = 1.0;
        let n = 32;
        let evenly_spaced_points = evenly_spaced_ellipse_points(a, b, n);
        let naively_spaced_points = naively_spaced_ellipse_points(a, b, n);
        assert!(distance_range(&naively_spaced_points) > distance_range(&evenly_spaced_points));
    }

    #[test]
    fn test_evenly_spaced_points_2() {
        let a = 1.0;
        let b = 1.5;
        let n = 32;
        let evenly_spaced_points = evenly_spaced_ellipse_points(a, b, n);
        let naively_spaced_points = naively_spaced_ellipse_points(a, b, n);
        assert!(distance_range(&naively_spaced_points) > distance_range(&evenly_spaced_points));
    }

    #[test]
    fn test_ellipse_circumference() {
        let width = 2.0;
        let height = 2.0;
        assert_abs_diff_eq!(ellipse_circumference(width, height), f64::consts::TAU * width);
    }
}