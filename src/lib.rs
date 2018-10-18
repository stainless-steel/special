//! [Special functions][1].
//!
//! [1]: https://en.wikipedia.org/wiki/Special_functions

#[cfg(test)]
extern crate assert;

extern crate libc;

mod bessel;
mod beta;
mod error;
mod gamma;
mod math;

pub use bessel::Bessel;
pub use beta::Beta;
pub use error::Error;
pub use gamma::Gamma;
