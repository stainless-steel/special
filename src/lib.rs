//! [Special functions][1].
//!
//! [1]: https://en.wikipedia.org/wiki/Special_functions

#[cfg(test)]
extern crate assert;

extern crate libc;

mod beta;
mod gamma;

#[link_name = "m"]
mod m;

pub use beta::{inc_beta, inv_inc_beta, ln_beta};
pub use gamma::{digamma, ln_gamma};

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
