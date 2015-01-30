//! [Special functions][1].
//!
//! [1]: https://en.wikipedia.org/wiki/Special_functions

#![feature(libc)]

#[cfg(test)]
#[macro_use]
extern crate assert;

#[cfg(test)]
extern crate test;

extern crate libc;

pub use beta::{ln_beta, inc_beta, inv_inc_beta};

mod beta;

#[link_name = "m"]
mod m {
    use libc::{c_double, c_int};

    extern {
        pub fn erf(x: c_double) -> c_double;
        pub fn erfc(x: c_double) -> c_double;
        pub fn lgamma_r(x: c_double, sign: &mut c_int) -> c_double;
    }
}

/// Compute the error function.
#[inline]
pub fn erf(x: f64) -> f64 {
    unsafe { m::erf(x) }
}

/// Compute the complementary error function.
#[inline]
pub fn erfc(x: f64) -> f64 {
    unsafe { m::erfc(x) }
}

/// Compute the natural logarithm of the gamma function.
#[inline]
pub fn ln_gamma(x: f64) -> (f64, i32) {
    let mut sign: i32 = 0;
    let value = unsafe { m::lgamma_r(x, &mut sign) };
    (value, sign)
}
