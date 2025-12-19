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
mod elliptic;
mod error;
mod gamma;
mod lambert_w;
mod primitive;

pub use crate::lambert_w::LambertW;
pub use beta::Beta;
pub use elliptic::Elliptic;
pub use error::Error;
pub use gamma::Gamma;
pub use primitive::Primitive;
