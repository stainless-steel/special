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
                unsafe { m::erf(self as f64) as Self }
            }

            #[inline]
            fn erfc(self) -> Self {
                unsafe { m::erfc(self as f64) as Self }
            }

            #[inline]
            fn gamma(self) -> Self {
                unsafe { m::tgamma(self as f64) as Self }
            }

            #[inline]
            fn ln_gamma(self) -> (Self, i32) {
                let mut sign: i32 = 0;
                let value = unsafe { m::lgamma(self as f64, &mut sign) as Self };
                (value, sign)
            }
        }
    )*);
}

implement!(f32, f64);

#[cfg(test)]
mod tests {
    use assert;

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

    #[test]
    fn gamma() {
        use std::f64::INFINITY;

        let x = vec![0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0];
        let y = vec![
                         INFINITY, 1.772453850905516e+00, 1.000000000000000e+00,
            8.862269254527581e-01, 1.000000000000000e+00, 1.329340388179137e+00,
            2.000000000000000e+00, 3.323350970447843e+00, 6.000000000000000e+00,
            1.163172839656745e+01, 2.400000000000000e+01, 5.234277778455353e+01,
            1.200000000000000e+02, 2.878852778150444e+02, 7.200000000000000e+02,
        ];

        let z = x.iter().map(|&x| x.gamma()).collect::<Vec<_>>();
        assert::close(&z, &y, 1e-12);
    }
}
