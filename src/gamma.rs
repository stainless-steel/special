use m;

/// Compute the real-valued digamma function.
///
/// The formula is as follows:
///
/// ```math
///        d ln(Γ(x))
/// ψ(x) = ----------
///            dx
/// ```
///
/// where Γ is the gamma function. The computation is based on an approximation
/// as described in the reference below.
///
/// ## Examples
///
/// ```
/// use special::digamma;
///
/// const EULER_MASCHERONI: f64 = 0.57721566490153286060651209008240243104215933593992;
/// assert!((digamma(1.0) + EULER_MASCHERONI).abs() < 1e-15);
/// ```
///
/// ## References
///
/// 1. M. J. Beal, Variational algorithms for approximate Bayesian inference.
///    University of London, 2003, pp. 265–266.
pub fn digamma(x: f64)-> f64 {
    macro_rules! evaluate_polynomial(
        ($x:expr, $coefficients:expr) => (
            $coefficients.iter().rev().fold(0.0, |sum, &c| $x * sum + c)
        );
    );

    if x <= 8.0 {
        return digamma(x + 1.0) - x.recip();
    }

    let inv_x = x.recip();
    let inv_x_2 = inv_x * inv_x;
    x.ln() - 0.5 * inv_x - inv_x_2 * evaluate_polynomial!(inv_x_2, [
        1.0 / 12.0, -1.0 / 120.0, 1.0 / 252.0, -1.0 / 240.0,
        5.0 / 660.0, -691.0 / 32760.0, 1.0 / 12.0, -3617.0 / 8160.0,
    ])
}

/// Compute the regularized lower incomplete gamma function.
///
/// The formula is as follows:
///
/// ```math
///           γ(x, p)    1   x
/// P(x, p) = ------- = ---- ∫ t^(p-1) e^(-t) dt
///            Γ(p)     Γ(p) 0
/// ```
///
/// where γ is the incomplete lower gamma function, and Γ is the complete gamma
/// function.
///
/// The code is based on a [C implementation][1] by John Burkardt. The original
/// algorithm was published in Applied Statistics and is known as
/// [Algorithm AS 239][2].
///
/// [1]: http://people.sc.fsu.edu/~jburkardt/c_src/asa239/asa239.html
/// [2]: http://www.jstor.org/stable/2347328
pub fn inc_gamma(x: f64, p: f64) -> f64 {
    debug_assert!(x >= 0.0 && p > 0.0);

    const ELIMIT: f64 = -88.0;
    const OFLO: f64 = 1.0e+37;
    const TOL: f64 = 1.0e-14;
    const XBIG: f64 = 1.0e+08;

    if x == 0.0 {
        return 0.0;
    }

    // For `p ≥ 1000`, the original algorithm uses an approximation shown below.
    // However, it introduces a substantial accuracy loss.
    //
    // ```
    // use std::f64::consts::FRAC_1_SQRT_2;
    //
    // const PLIMIT: f64 = 1000.0;
    //
    // if PLIMIT < p {
    //     let pn1 = 3.0 * p.sqrt() * ((x / p).powf(1.0 / 3.0) + 1.0 / (9.0 * p) - 1.0);
    //     return 0.5 * (1.0 + unsafe { m::erf(FRAC_1_SQRT_2 * pn1) });
    // }
    // ```

    if XBIG < x {
        return 1.0;
    }

    if x <= 1.0 || x < p {
        let mut arg = p * x.ln() - x - ln_gamma(p + 1.0).0;
        let mut c = 1.0;
        let mut value = 1.0;
        let mut a = p;

        loop {
            a += 1.0;
            c *= x / a;
            value += c;

            if c <= TOL {
                break;
            }
        }

        arg += value.ln();

        if ELIMIT <= arg {
            return arg.exp();
        } else {
            return 0.0;
        };
    } else {
        let mut arg = p * x.ln() - x - ln_gamma(p).0;
        let mut a = 1.0 - p;
        let mut b = a + x + 1.0;
        let mut c = 0.0;
        let mut pn1 = 1.0;
        let mut pn2 = x;
        let mut pn3 = x + 1.0;
        let mut pn4 = x * b;
        let mut value = pn3 / pn4;

        loop {
            a += 1.0;
            b += 2.0;
            c += 1.0;
            let an = a * c;
            let pn5 = b * pn3 - an * pn1;
            let pn6 = b * pn4 - an * pn2;

            if pn6 != 0.0 {
                let rn = pn5 / pn6;
                if (value - rn).abs() <= TOL.min(TOL * rn) {
                    break;
                }
                value = rn;
            }

            pn1 = pn3;
            pn2 = pn4;
            pn3 = pn5;
            pn4 = pn6;

            if OFLO <= pn5.abs() {
                pn1 /= OFLO;
                pn2 /= OFLO;
                pn3 /= OFLO;
                pn4 /= OFLO;
            }
        }

        arg += value.ln();

        if ELIMIT <= arg {
            return 1.0 - arg.exp();
        } else {
            return 1.0;
        }
    }
}

/// Compute the natural logarithm of the gamma function.
#[inline]
pub fn ln_gamma(x: f64) -> (f64, i32) {
    let mut sign: i32 = 0;
    let value = unsafe { m::lgamma(x, &mut sign) };
    (value, sign)
}

#[cfg(test)]
mod tests {
    use assert;

    #[test]
    fn digamma() {
        use std::f64::consts::{FRAC_PI_2, LN_2};
        const EULER_MASCHERONI: f64 = 0.57721566490153286060651209008240243104215933593992;
        assert_eq!(-FRAC_PI_2 - 3.0 * LN_2 - EULER_MASCHERONI, super::digamma(0.25));
    }

    #[test]
    fn inc_gamma_small_p() {
        let p = 4.2;
        let x = vec![
            0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5,
            5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5,
        ];
        let y = vec![
            0.000000000000000e+00, 1.118818097571853e-03, 1.386936454406691e-02,
            5.173931341260211e-02, 1.186345016135447e-01, 2.092219100370761e-01,
            3.137074857845146e-01, 4.219736975484903e-01, 5.258832878449179e-01,
            6.200417684429613e-01, 7.016346974530381e-01, 7.698559854309833e-01,
            8.252531208548146e-01, 8.691540627528154e-01, 9.032348314835232e-01,
            9.292289044193606e-01, 9.487539837941871e-01, 9.632249236076943e-01,
            9.738240752336722e-01, 9.815062931491536e-01,
        ];

        let z = x.iter().map(|&x| super::inc_gamma(x, p)).collect::<Vec<_>>();
        assert::close(&z, &y, 1e-14);
    }

    #[test]
    fn inc_gamma_large_p() {
        let p = 1500.0;
        let x = vec![
            1400.0, 1410.0, 1420.0, 1430.0, 1440.0, 1450.0, 1460.0, 1470.0, 1480.0, 1490.0,
            1500.0, 1510.0, 1520.0, 1530.0, 1540.0, 1550.0, 1560.0, 1570.0, 1580.0, 1590.0,
        ];
        let y = vec![
            4.231080348517120e-03, 9.056183893278287e-03, 1.809200095269094e-02,
            3.380047876608895e-02, 5.917825204288418e-02, 9.731671903009434e-02,
            1.506860564044825e-01, 2.202940362119810e-01, 3.049926190158564e-01,
            4.012305141433336e-01, 5.034335611603484e-01, 6.049687028759408e-01,
            6.994146092246781e-01, 7.817404371226420e-01, 8.490443067241588e-01,
            9.006922020203059e-01, 9.379249194843459e-01, 9.631597573132408e-01,
            9.792520392189441e-01, 9.889149370075309e-01,
        ];

        let z = x.iter().map(|&x| super::inc_gamma(x, p)).collect::<Vec<_>>();
        assert::close(&z, &y, 1e-12);
    }
}
