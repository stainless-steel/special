//! The library provides [special functions](
//! https://en.wikipedia.org/wiki/Special_functions).

#![feature(macro_rules)]

extern crate libc;

mod cmath {
    use libc::{c_double, c_int};

    #[link_name = "m"]
    extern {
        pub fn erf(x: c_double) -> c_double;
        pub fn erfc(x: c_double) -> c_double;
        pub fn exp(x: c_double) -> c_double;
        pub fn lgamma_r(x: c_double, sign: &mut c_int) -> c_double;
    }
}

/// Computes the error function.
pub fn erf(x: f64) -> f64 {
    unsafe { cmath::erf(x) }
}

/// Computes the complementary error function.
pub fn erfc(x: f64) -> f64 {
    unsafe { cmath::erfc(x) }
}

/// Computes the exponential function.
pub fn exp(x: f64) -> f64 {
    unsafe { cmath::exp(x) }
}

/// Computes the natural logarithm of the beta function.
pub fn lbeta(x: f64, y: f64) -> f64 {
    lgamma(x) + lgamma(y) - lgamma(x + y)
}

/// Computes the natural logarithm of the gamma function.
pub fn lgamma(x: f64) -> f64 {
    unsafe {
        let mut sign: libc::c_int = 0;
        cmath::lgamma_r(x, &mut sign)
    }
}

#[cfg(test)]
mod tests {
    macro_rules! assert_close(
        ($x:expr, $y:expr) => ({
            let eps: f64 = ::std::f64::EPSILON.sqrt();
            let x: Vec<f64> = $x;
            let y: Vec<f64> = $y;
            for i in range(0u, x.len()) {
                assert!(::std::num::abs(x[i] - y[i]) < eps,
                       "expected {:e} ~ {:e}", x[i], y[i]);
            }
        });
    )

    macro_rules! apply(
        ($x:expr, $func:expr) => (
            $x.iter().map(|&x| $func(x)).collect()
        );
    )

    macro_rules! apply2(
        ($x:expr, $y:expr, $func:expr) => (
            $x.iter().zip($y.iter()).map(|(&x, &y)| $func(x, y)).collect()
        );
    )

    #[test]
    fn erf() {
        let x = vec![-3.0, -1.0, 1.0, 3.0];
        let y = vec![-0.9999779095030014, -0.8427007929497149,
                     0.8427007929497149, 0.9999779095030014];
        assert_close!(apply!(x, super::erf), y);
    }

    #[test]
    fn erfc() {
        let x = vec![-3.0, -1.0, 1.0, 3.0];
        let y = vec![1.9999779095030015, 1.8427007929497150,
                     0.1572992070502851, 0.0000220904969986];
        assert_close!(apply!(x, super::erfc), y);
    }

    #[test]
    fn exp() {
        let x = vec![-2.0, -1.0, 1.0, 2.0];
        let y = vec![0.1353352832366127, 0.3678794411714423,
                     2.7182818284590455, 7.3890560989306504];
        assert_close!(apply!(x, super::exp), y);
    }

    #[test]
    fn lbeta() {
        let x = vec![0.25, 0.50, 0.75, 1.00];
        let y = vec![0.50, 0.75, 1.00, 1.25];
        let z = vec![1.6571065161914820, 0.8739177307778084,
                     0.2876820724517809, -0.2231435513142098];
        assert_close!(apply2!(x, y, super::lbeta), z);
    }

    #[test]
    fn lgamma() {
        let x = vec![0.25, 0.50, 0.75, 1.00];
        let y = vec![1.2880225246980774, 0.5723649429247000,
                     0.2032809514312953, 0.0000000000000000];
        assert_close!(apply!(x, super::lgamma), y);
    }
}
