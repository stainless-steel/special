//! [Special functions][1].
//!
//! [1]: https://en.wikipedia.org/wiki/Special_functions

#[cfg(test)]
extern crate assert;


extern crate libc;
extern crate num;

pub use beta::{ln_beta, inc_beta, inv_inc_beta};

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

// Use Horners method to evaluate polynomials at `x0`
macro_rules! eval_poly {
    ( $x0:expr; $( $p:expr ),* ) => {
        {
            let mut coefs = vec![$($p,)*];
            coefs.reverse();
            let mut sum = 0.0;
            for c in coefs.iter() {
                sum = sum * $x0 + c;
            }
            sum
        }
    }
}

/// Compute the real valued digamma function. `d/dx ln Gamma(x)`.
/// Appriximation based on Beal, Matthew J. (2003). *Variational
/// Algorithms for Approximate Bayesian Inference* (PhD thesis). The
/// Gatsby Computational Neuroscience Unit, University College
/// London. pp. 265–266
///
/// # Examples
///
/// ```
/// use special::digamma;
///
/// let euler_mascheroni = 0.5772156649015325;
/// assert_eq!(-euler_mascheroni, digamma(1.0));
/// ```
#[inline]
pub fn digamma(x: f64)-> f64 {
    if x > 8.0 {
        let inv_x = x.recip();
        let inv_x_e2 = inv_x.powi(2);
        x.ln() - 0.5*inv_x - inv_x_e2 * eval_poly!(
            inv_x_e2;
            // Coefficients are `Bernoulli[2n] / 2n`
            1.0/12.0, -1.0/120.0, 1.0/252.0, -1.0/240.0,
            5.0/660.0, -691.0/32760.0, 1.0/12.0, -3617.0/8160.0
        )
    } else {
        // Shift by recurrence
        digamma(x + 1.0) - x.recip()
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_digamma() {
        use std::f64::consts::{ FRAC_PI_2, LN_2 };
        let euler_mascheroni = 0.577215664901533;

        assert_eq!(-FRAC_PI_2 - 3.0 * LN_2 - euler_mascheroni, digamma(0.25));
    }
}
