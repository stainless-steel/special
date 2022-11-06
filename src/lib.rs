//! [Special functions][1].
//!
//! [1]: https://en.wikipedia.org/wiki/Special_functions

#![no_std]

#[cfg(test)]
extern crate assert;

extern crate alloc;
extern crate libm;

mod basic;
mod beta;
mod error;
mod gamma;

pub use basic::Basic;
pub use beta::Beta;
pub use error::Error;
pub use gamma::Gamma;
