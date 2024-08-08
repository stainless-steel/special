//! [Special functions][1].
//!
//! [1]: https://en.wikipedia.org/wiki/Special_functions

#![cfg_attr(feature = "no_std", no_std)]
#![allow(clippy::excessive_precision)]

#[cfg(test)]
extern crate assert;

#[cfg(test)]
extern crate alloc;

mod beta;
mod error;
mod gamma;
mod omega;
mod primitive;

pub use beta::Beta;
pub use error::Error;
pub use gamma::Gamma;
pub use omega::Omega;
pub use primitive::Primitive;
