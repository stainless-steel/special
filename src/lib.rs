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
#[cfg(feature = "elliptic")]
mod elliptic;
mod error;
mod gamma;
#[cfg(feature = "lambert_w")]
mod lambert_w;
mod primitive;

pub use crate::beta::Beta;
#[cfg(feature = "elliptic")]
pub use crate::elliptic::Elliptic;
pub use crate::error::Error;
pub use crate::gamma::Gamma;
#[cfg(feature = "lambert_w")]
pub use crate::lambert_w::LambertW;
pub use crate::primitive::Primitive;
