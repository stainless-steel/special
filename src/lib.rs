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

pub use beta::{inc_beta, inv_inc_beta, ln_beta};
pub use error::Error;
pub use gamma::{Gamma, inc_gamma};

/// Special functions.
pub trait Special where Self: Error + Gamma {
}

impl<T> Special for T where T: Error + Gamma {
}
