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

/// Computes `x` raised to the exponent `y`.
pub fn pow(x: f64, y: f64) -> f64 {
    unsafe { m::pow(x, y) }
}

/// Computes the positive square root of the argument.
pub fn sqrt(x: f64) -> f64 {
    unsafe { m::sqrt(x) }
}
