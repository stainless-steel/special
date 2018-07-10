use math;
use std::f64;

const SQRT_PI: f64 = 1.772453850905515881919427556567825376987457275391f64;

/// Error functions.
pub trait Error {
    /// Compute the error function.
    fn erf(self) -> Self;

    /// Compute the complementary error function.
    fn erfc(self) -> Self;

    /// Compute the inverse error function.
    fn erfinv(self) -> Self;
}

// Code translated from C:
// https://scistatcalc.blogspot.com/2013/09/numerical-estimate-of-inverse-error.html
fn erfinv(z: f64) -> f64 {
    let mut w: f64 = -((1.0 - z) * (1.0 + z)).ln();
    let mut p: f64;

    if w < 5.0 {
        w -= 2.5;
        p = 2.81022636e-08;
        p = 3.43273939e-07 + p * w;
        p = -3.5233877e-06 + p * w;
        p = -4.39150654e-06 + p * w;
        p = 0.00021858087 + p * w;
        p = -0.00125372503 + p * w;
        p = -0.00417768164 + p * w;
        p = 0.246640727 + p * w;
        p = 1.50140941 + p * w;
    } else {
        w = w.sqrt() - 3.0;
        p = -0.000200214257;
        p = 0.000100950558 + p * w;
        p = 0.00134934322 + p * w;
        p = -0.00367342844 + p * w;
        p = 0.00573950773 + p * w;
        p = -0.0076224613 + p * w;
        p = 0.00943887047 + p * w;
        p = 1.00167406 + p * w;
        p = 2.83297682 + p * w;
    }

    let res_ra = p * z; // assign to rational estimate variable

    // Halley's method to refine estimate of inverse erf
    let fx = unsafe { math::erf(res_ra) } - z;
    let df = 2.0 / SQRT_PI * (-(res_ra * res_ra)).exp();
    let d2f = -2.0 * res_ra * df;

    res_ra - (2.0 * fx * df) / ((2.0 * df * df) - (fx * d2f))
}

macro_rules! implement {
    ($kind:ty) => {
        impl Error for $kind {
            #[inline]
            fn erf(self) -> Self {
                unsafe { math::erf(self as f64) as Self }
            }

            #[inline]
            fn erfc(self) -> Self {
                unsafe { math::erfc(self as f64) as Self }
            }

            #[inline]
            fn erfinv(self) -> Self {
                erfinv(self as f64) as Self
            }
        }
    };
}

implement!(f32);
implement!(f64);

#[cfg(test)]
mod tests {
    use super::*;
    use assert;

    // Test values calculated from scipy
    #[test]
    fn erfinv_of_zero_should_be_zero() {
        assert::close(erfinv(0.0), 0.0, 1E-12);
    }

    #[test]
    fn erfinv_positive_value_mid_domain() {
        assert::close(erfinv(0.5), 0.47693627620446982, 1E-12);
    }

    #[test]
    fn erfinv_positive_value_near_center() {
        assert::close(erfinv(0.121), 0.10764782605515244, 1E-12);
    }

    #[test]
    fn erfinv_negative_value_near_lower_bound() {
        assert::close(erfinv(-0.99), -1.8213863677184492, 1E-12);
    }

    #[test]
    fn erfinv_negative_value_nearer_lower_bound() {
        assert::close(erfinv(-0.999), -2.3267537655135242, 1E-12);
    }
}
