//! [Special functions][1].
//!
//! [1]: https://en.wikipedia.org/wiki/Special_functions

#![no_std]

#[cfg(test)]
extern crate assert;

extern crate alloc;
extern crate libm;

mod beta;
mod error;
mod gamma;
mod primitive;

pub use beta::Beta;
pub use error::Error;
pub use gamma::Gamma;
pub use primitive::Primitive;
