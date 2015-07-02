//! [Special functions][1].
//!
//! [1]: https://en.wikipedia.org/wiki/Special_functions

#[cfg(test)]
extern crate assert;

extern crate libc;

mod beta;

#[link_name = "m"]
mod m {
    use libc::{c_double, c_int};

    extern {
        pub fn erf(x: c_double) -> c_double;
        pub fn erfc(x: c_double) -> c_double;
        pub fn lgamma_r(x: c_double, sign: &mut c_int) -> c_double;
    }
}

pub use beta::{ln_beta, inc_beta, inv_inc_beta};

/// Compute the real-valued digamma function.
///
/// The digamma function is defined by
///
/// ```math
///        d ln(Γ(x))
/// ψ(x) = ----------
///            dx
/// ```
///
/// where Γ is the Gamma function. The computation is based on an approximation
/// as described in the reference below.
///
/// ## Examples
///
/// ```
/// use special::digamma;
///
/// let euler_mascheroni = 0.5772156649015325;
/// assert_eq!(-euler_mascheroni, digamma(1.0));
/// ```
///
/// ## References
///
/// 1. M. J. Beal, Variational algorithms for approximate Bayesian inference.
///    University of London, 2003, pp. 265–266.
pub fn digamma(x: f64)-> f64 {
    macro_rules! evaluate_polynom(
        ($x:expr, $coefficients:expr) => (
            $coefficients.iter().rev().fold(0.0, |sum, &c| $x * sum + c)
        );
    );

    if x <= 8.0 {
        return digamma(x + 1.0) - x.recip();
    }

    let inv_x = x.recip();
    let inv_x_2 = inv_x * inv_x;
    x.ln() - 0.5 * inv_x - inv_x_2 * evaluate_polynom!(inv_x_2, [
        1.0 / 12.0, -1.0 / 120.0, 1.0 / 252.0, -1.0 / 240.0,
        5.0 / 660.0, -691.0 / 32760.0, 1.0 / 12.0, -3617.0 / 8160.0,
    ])
}

/// Compute the error function.
#[inline]
pub fn erf(x: f64) -> f64 {
    unsafe { m::erf(x) }
}

/// Compute the complementary error function.
#[inline]
pub fn erfc(x: f64) -> f64 {
    unsafe { m::erfc(x) }
}

/// Compute the natural logarithm of the gamma function.
#[inline]
pub fn ln_gamma(x: f64) -> (f64, i32) {
    let mut sign: i32 = 0;
    let value = unsafe { m::lgamma_r(x, &mut sign) };
    (value, sign)
}

#[cfg(test)]
mod tests {
    #[test]
    fn digamma() {
        use std::f64::consts::{FRAC_PI_2, LN_2};
        let euler_mascheroni = 0.577215664901533;
        assert_eq!(-FRAC_PI_2 - 3.0 * LN_2 - euler_mascheroni, super::digamma(0.25));
    }
}
