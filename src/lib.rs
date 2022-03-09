//! [Special functions][1].
//!
//! [1]: https://en.wikipedia.org/wiki/Special_functions

#![cfg_attr(not(feature = "std"), no_std)]

#[cfg(feature = "std")]
extern crate core;

extern crate alloc;

#[cfg(test)]
extern crate assert;

#[cfg(feature = "system_math")]
extern crate libc;

#[cfg(not(feature = "system_math"))]
extern crate libm;

mod beta;
mod consts;
mod error;
mod gamma;
mod math;

pub use beta::Beta;
pub use error::Error;
pub use gamma::Gamma;
pub use math::FloatExt;
