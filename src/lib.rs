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
pub use gamma::{digamma, gamma, inc_gamma, ln_gamma};

/// Special functions.
pub trait Special {
    /// Compute the error function.
    fn erf(self) -> Self;

    /// Compute the complementary error function.
    fn erfc(self) -> Self;
}

macro_rules! implement {
    ($($kind:ty),*) => ($(
        impl Special for $kind {
            #[inline]
            fn erf(self) -> $kind {
                unsafe { m::erf(self as f64) as $kind }
            }

            #[inline]
            fn erfc(self) -> $kind {
                unsafe { m::erfc(self as f64) as $kind }
            }
        }
    )*);
}

implement!(f32, f64);

#[cfg(test)]
mod tests {
    use Special;

    #[test]
    pub fn erf() {
        0f32.erf();
        0f64.erf();
    }

    #[test]
    pub fn erfc() {
        0f32.erfc();
        0f64.erfc();
    }
}
