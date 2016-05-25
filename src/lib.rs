//! [Special functions][1].
//!
//! [1]: https://en.wikipedia.org/wiki/Special_functions

#[cfg(test)]
extern crate assert;

extern crate libc;

mod beta;
mod gamma;
mod math;

pub use beta::{inc_beta, inv_inc_beta, ln_beta};
pub use gamma::{digamma, inc_gamma};

/// Special functions.
pub trait Special where Self: Sized {
    /// Compute the error function.
    fn erf(self) -> Self;

    /// Compute the complementary error function.
    fn erfc(self) -> Self;

    /// Compute the gamma function.
    fn gamma(self) -> Self;

    /// Compute the natural logarithm of the gamma function.
    fn ln_gamma(self) -> (Self, i32);
}

macro_rules! implement {
    ($($kind:ty),*) => ($(
        impl Special for $kind {
            #[inline]
            fn erf(self) -> Self {
                unsafe { math::erf(self as f64) as Self }
            }

            #[inline]
            fn erfc(self) -> Self {
                unsafe { math::erfc(self as f64) as Self }
            }

            #[inline]
            fn gamma(self) -> Self {
                unsafe { math::tgamma(self as f64) as Self }
            }

            #[inline]
            fn ln_gamma(self) -> (Self, i32) {
                let mut sign: i32 = 0;
                let value = unsafe { math::lgamma(self as f64, &mut sign) as Self };
                (value, sign)
            }
        }
    )*);
}

implement!(f32, f64);
