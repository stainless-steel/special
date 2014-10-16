//! The library provides [special functions][1].
//!
//! [1]: https://en.wikipedia.org/wiki/Special_functions

#![feature(macro_rules)]

extern crate libc;

pub use beta::{log_beta, inc_beta, inv_inc_beta};

mod beta;
mod m;

/// Computes the error function.
pub fn erf(x: f64) -> f64 {
    unsafe { m::erf(x) }
}

/// Computes the complementary error function.
pub fn erfc(x: f64) -> f64 {
    unsafe { m::erfc(x) }
}

/// Computes the exponential function.
pub fn exp(x: f64) -> f64 {
    unsafe { m::exp(x) }
}

/// Computes the natural logarithm of the gamma function.
pub fn log_gamma(x: f64) -> f64 {
    unsafe {
        let mut sign: i32 = 0;
        m::lgamma_r(x, &mut sign)
    }
}

/// Computes the natural logarithm.
pub fn log(x: f64) -> f64 {
    unsafe { m::log(x) }
}

#[cfg(test)]
mod tests {
    macro_rules! assert_close(
        ($x:expr, $y:expr) => ({
            let eps: f64 = ::std::f64::EPSILON.sqrt();
            let x: Vec<f64> = $x;
            let y: Vec<f64> = $y;
            for i in range(0u, x.len()) {
                assert!(::std::num::abs(x[i] - y[i]) < eps,
                       "expected {:e} ~ {:e}", x[i], y[i]);
            }
        });
    )

    macro_rules! apply(
        ($x:expr, $func:expr) => (
            $x.iter().map(|&x| $func(x)).collect()
        );
    )

    #[test]
    fn erf() {
        let x = vec![-3.0, -1.0, 1.0, 3.0];
        let y = vec![-0.9999779095030014, -0.8427007929497149,
                     0.8427007929497149, 0.9999779095030014];
        assert_close!(apply!(x, super::erf), y);
    }

    #[test]
    fn erfc() {
        let x = vec![-3.0, -1.0, 1.0, 3.0];
        let y = vec![1.9999779095030015, 1.8427007929497150,
                     0.1572992070502851, 0.0000220904969986];
        assert_close!(apply!(x, super::erfc), y);
    }

    #[test]
    fn exp() {
        let x = vec![-2.0, -1.0, 1.0, 2.0];
        let y = vec![0.1353352832366127, 0.3678794411714423,
                     2.7182818284590455, 7.3890560989306504];
        assert_close!(apply!(x, super::exp), y);
    }

    #[test]
    fn log_gamma() {
        let x = vec![0.25, 0.50, 0.75, 1.00];
        let y = vec![1.2880225246980774, 0.5723649429247000,
                     0.2032809514312953, 0.0000000000000000];
        assert_close!(apply!(x, super::log_gamma), y);
    }

    #[test]
    fn log() {
        let x = vec![0.5, 1.0, 1.5, 2.0];
        let y = vec![-0.6931471805599453, 0.0000000000000000,
                     0.4054651081081644, 0.6931471805599453];
        assert_close!(apply!(x, super::log), y);
    }
}
