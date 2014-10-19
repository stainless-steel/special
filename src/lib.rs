//! The library provides [special functions][1].
//!
//! [1]: https://en.wikipedia.org/wiki/Special_functions

#![feature(macro_rules)]

extern crate libc;

pub use beta::{log_beta, inc_beta, inv_inc_beta};

mod beta;

#[link_name = "m"]
mod m {
    use libc::{c_double, c_int};

    extern {
        pub fn cos(x: c_double) -> c_double;
        pub fn erf(x: c_double) -> c_double;
        pub fn erfc(x: c_double) -> c_double;
        pub fn exp(x: c_double) -> c_double;
        pub fn tgamma(x: c_double) -> c_double;
        pub fn lgamma_r(x: c_double, sign: &mut c_int) -> c_double;
        pub fn log(x: c_double) -> c_double;
        pub fn pow(x: c_double, y: c_double) -> c_double;
        pub fn sin(x: c_double) -> c_double;
        pub fn sqrt(x: c_double) -> c_double;
    }
}

/// Computes the cosine.
#[inline]
pub fn cos(x: f64) -> f64 {
    unsafe { m::cos(x) }
}

/// Computes the error function.
#[inline]
pub fn erf(x: f64) -> f64 {
    unsafe { m::erf(x) }
}

/// Computes the complementary error function.
#[inline]
pub fn erfc(x: f64) -> f64 {
    unsafe { m::erfc(x) }
}

/// Computes the exponential function.
#[inline]
pub fn exp(x: f64) -> f64 {
    unsafe { m::exp(x) }
}

/// Computes the gamma function.
#[inline]
pub fn gamma(x: f64) -> f64 {
    unsafe { m::tgamma(x) }
}

/// Computes the natural logarithm of the gamma function.
#[inline]
pub fn log_gamma(x: f64) -> (f64, i32) {
    let mut sign: i32 = 0;
    let value = unsafe { m::lgamma_r(x, &mut sign) };
    (value, sign)
}

/// Computes the natural logarithm.
///
/// The domain is positive real numbers.
#[inline]
pub fn log(x: f64) -> f64 {
    unsafe { m::log(x) }
}

/// Computes a base raised to an exponent.
#[inline]
pub fn pow(base: f64, exp: f64) -> f64 {
    unsafe { m::pow(base, exp) }
}

/// Computes the sine.
#[inline]
pub fn sin(x: f64) -> f64 {
    unsafe { m::sin(x) }
}

/// Computes the square root function.
///
/// The domain is non-negative real numbers.
#[inline]
pub fn sqrt(x: f64) -> f64 {
    unsafe { m::sqrt(x) }
}

#[cfg(test)]
mod test {
    #[test]
    fn gamma() {
        use super::{exp, gamma, log_gamma};

        let x = -2.5;

        let gamma1 = -0.9453087204829418;
        let gamma2 = gamma(x);
        let gamma3 = {
            let (value, sign) = log_gamma(x);
            (sign as f64) * exp(value)
        };

        assert_eq!(gamma1, gamma2);
        assert_eq!(gamma1, gamma3);
    }

}
