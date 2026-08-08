use crate::primitive::Primitive;

// The value is found empirically. It must be even; otherwise, it lands on the wrong phase of a
// 1-ulp oscillation at the deepest tail points.
const HALLEY_ITERATIONS: usize = 8;

/// Error functions.
pub trait Error {
    /// Compute the error function.
    fn error(self) -> Self;

    /// Compute the complementary error function.
    fn compl_error(self) -> Self;

    /// Compute the inverse of the error function.
    ///
    /// The implementation is based on a [C implementation][1] by Alijah Ahmed.
    ///
    /// [1]: https://scistatcalc.blogspot.com/2013/09/numerical-estimate-of-inverse-error.html
    fn inv_error(self) -> Self;
}

macro_rules! implement {
    ($kind:ty) => {
        impl Error for $kind {
            #[inline]
            fn error(self) -> Self {
                Primitive::erf(self)
            }

            #[inline]
            fn compl_error(self) -> Self {
                Primitive::erfc(self)
            }

            fn inv_error(self) -> Self {
                const SQRT_PI: $kind = 1.772453850905515881919427556567825376987457275391;

                let mut w: $kind = -((1.0 - self) * (1.0 + self)).ln();
                let mut p: $kind;

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

                    let res_ra = p * self;
                    let fx: Self = res_ra.error() - self;
                    let df = 2.0 / SQRT_PI as $kind * (-(res_ra * res_ra)).exp();
                    let d2f = -2.0 * res_ra * df;

                    res_ra - (2.0 * fx * df) / ((2.0 * df * df) - (fx * d2f))
                } else if self == 1.0 {
                    return <$kind>::INFINITY;
                } else if self == -1.0 {
                    return <$kind>::NEG_INFINITY;
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

                    let mut res_ra = p * self;
                    for _ in 0..HALLEY_ITERATIONS {
                        let fx: Self = if res_ra >= 0.0 {
                            (1.0 - self) - res_ra.compl_error()
                        } else {
                            (-res_ra).compl_error() - (1.0 + self)
                        };
                        let df = 2.0 / SQRT_PI as $kind * (-(res_ra * res_ra)).exp();
                        let d2f = -2.0 * res_ra * df;
                        let next = res_ra - (2.0 * fx * df) / ((2.0 * df * df) - (fx * d2f));
                        if next == res_ra {
                            break;
                        }
                        res_ra = next;
                    }
                    res_ra
                }
            }
        }
    };
}

implement!(f32);
implement!(f64);

#[cfg(test)]
mod tests {
    use assert;

    use super::*;

    #[test]
    fn inv_error_negative() {
        assert::close(-0.99.inv_error(), -1.8213863677184492, 1e-12);
        assert::close(-0.999.inv_error(), -2.3267537655135242, 1e-12);
    }

    #[test]
    fn inv_error_positive() {
        assert::close(0.5.inv_error(), 0.47693627620446982, 1e-12);
        assert::close(0.121.inv_error(), 0.10764782605515244, 1e-12);
    }

    #[test]
    fn inv_error_zero() {
        assert::close(0.0.inv_error(), 0.0, 1e-12);
    }

    #[test]
    fn inv_error_deep_tail() {
        // Correctly rounded values of erfinv(1 - 2^-k) for the tail region
        // (mpmath at 65 and 90 digits, agreeing bitwise after rounding).
        const CASES: [(i32, u64); 19] = [
            (9, 0x4001855321d8668a),
            (10, 0x4002a6d8937b12d6),
            (11, 0x4003b9ddbc5d5000),
            (12, 0x4004c04ee987ec04),
            (13, 0x4005bbb5be25f0e8),
            (14, 0x4006ad52a4a63e03),
            (16, 0x40087726db4ac1d2),
            (18, 0x400a2441e0f076db),
            (20, 0x400bb95ab2ef976f),
            (24, 0x400ea8f95ae0aca8),
            (28, 0x4010ad1fec244a36),
            (32, 0x4011ed2bef15a0b5),
            (36, 0x401319301128afb8),
            (40, 0x4014347bf36fbae8),
            (44, 0x4015418dbe470875),
            (48, 0x40164253f064686f),
            (50, 0x4016be9874412264),
            (52, 0x40173856d153f081),
            (53, 0x4017744f8f74e94a),
        ];

        for &(k, bits) in &CASES {
            let y = 1.0 - (2.0_f64).powi(-k);
            let expected = f64::from_bits(bits);
            assert_eq!(y.inv_error(), expected);
            assert_eq!((-y).inv_error(), -expected);
        }

        assert_eq!(1.0_f64.inv_error(), f64::INFINITY);
        assert_eq!((-1.0_f64).inv_error(), f64::NEG_INFINITY);
    }

    #[test]
    fn inv_error_f32_tail() {
        const CASES: [(i32, u32); 5] = [
            (11, 0x401dceee),
            (13, 0x402dddae),
            (16, 0x4043b937),
            (21, 0x4063e060),
            (24, 0x407547cb),
        ];

        for &(k, bits) in &CASES {
            let y = 1.0_f32 - (2.0_f32).powi(-k);
            let expected = f32::from_bits(bits);
            assert_eq!(y.inv_error(), expected);
            assert_eq!((-y).inv_error(), -expected);
        }

        assert_eq!(1.0_f32.inv_error(), f32::INFINITY);
        assert_eq!((-1.0_f32).inv_error(), f32::NEG_INFINITY);
    }
}
