//! [Special functions][1].
//!
//! [1]: https://en.wikipedia.org/wiki/Special_functions

#[cfg(test)]
extern crate assert;

extern crate libc;

mod beta;
mod error;
mod gamma;
mod math;

pub use beta::Beta;
pub use error::Error;
pub use gamma::Gamma;

/// Special functions.
pub trait Special where Self: Beta + Error + Gamma {
}

impl<T> Special for T where T: Beta + Error + Gamma {
}
