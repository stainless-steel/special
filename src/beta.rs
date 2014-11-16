/// Compute the natural logarithm of the beta function.
///
/// The domain is `{(x, y): x > 0, y > 0}`.
pub fn ln_beta(x: f64, y: f64) -> f64 {
    use super::ln_gamma;

    let (a, _) = ln_gamma(x);
    let (b, _) = ln_gamma(y);
    let (c, _) = ln_gamma(x + y);

    a + b - c
}

/// Compute the incomplete beta function.
///
/// The code is based on a [C implementation][1] by John Burkardt. The
/// original algorithm was published in Applied Statistics and is known as
/// [Algorithm AS 63][2] and [Algorithm AS 109][3].
///
/// [1]: http://people.sc.fsu.edu/~jburkardt/c_src/asa109/asa109.html
/// [2]: http://www.jstor.org/stable/2346797
/// [3]: http://www.jstor.org/stable/2346887
pub fn inc_beta(x: f64, mut p: f64, mut q: f64, ln_beta: f64) -> f64 {
    // Algorithm AS 63
    // http://www.jstor.org/stable/2346797
    //
    // The function uses the method discussed by Soper (1921). If p is not less
    // than (p + q)x and the integral part of q + (1 - x)(p + q) is a positive
    // integer, say s, reductions are made up to s times “by parts” using the
    // recurrence relation
    //
    //                 Γ(p+q)
    // I(x, p, q) = ----------- x^p (1-x)^(q-1) + I(x, p+1, q-1)
    //              Γ(p+1) Γ(q)
    //
    // and then reductions are continued by “raising p” with the recurrence
    // relation
    //
    //                      Γ(p+q)
    // I(x, p+s, q-s) = --------------- x^(p+s) (1 - x)^(q-s) + I(x, p+s+1, q-s)
    //                  Γ(p+s+1) Γ(q-s)
    //
    // If s is not a positive integer, reductions are made only by “raising p.”
    // The process of reduction is terminated when the relative contribution to
    // the integral is not greater than the value of ACU. If p is less than
    // (p + q)x, I(1-x, q, p) is first calculated by the above procedure and
    // then I(x, p, q) is obtained from the relation
    //
    // I(x, p, q) = 1 - I(1-x, p, q).
    //
    // Soper (1921) demonstrated that the expansion of I(x, p, q) by “parts”
    // and “raising p” method as described above converges more rapidly than
    // any other series expansions.

    use std::num::Float;

    #[inline(always)]
    fn exp(x: f64) -> f64 { x.exp() }
    #[inline(always)]
    fn ln(x: f64) -> f64 { x.ln() }

    const ACU: f64 = 0.1e-14;

    if x <= 0.0 {
        return 0.0;
    }
    if 1.0 <= x {
        return 1.0;
    }

    let mut psq = p + q;

    let mut pbase;
    let mut qbase;

    let mut temp;

    let flip = p < psq * x;
    if flip {
        pbase = 1.0 - x;
        qbase = x;
        temp = q;
        q = p;
        p = temp;
    } else {
        pbase = x;
        qbase = 1.0 - x;
    }

    let mut term = 1.0;
    let mut ai = 1.0;

    let mut rx;
    let mut ns = (q + qbase * psq) as int;
    if ns == 0 {
        rx = pbase;
    } else {
        rx = pbase / qbase;
    }

    let mut alpha = 1.0;
    temp = q - ai;

    loop {
        term = term * temp * rx / (p + ai);

        alpha += term;

        temp = if term < 0.0 { -term } else { term };
        if temp <= ACU && temp <= ACU * alpha {
            break;
        }

        ai += 1.0;
        ns -= 1;

        if 0 < ns {
            temp = q - ai;
        } else if ns == 0 {
            temp = q - ai;
            rx = pbase;
        } else {
            temp = psq;
            psq += 1.0;
        }
    }

    // Remark AS R19 and Algorithm AS 109
    // http://www.jstor.org/stable/2346887
    alpha = alpha * exp(p * ln(pbase) + (q - 1.0) * ln(qbase) - ln_beta) / p;

    if flip { 1.0 - alpha } else { alpha }
}

/// Compute the inverse of the incomplete beta function.
///
/// The code is based on a [C implementation][1] by John Burkardt. The
/// original algorithm was published in Applied Statistics and is known as
/// [Algorithm AS 64][2] and [Algorithm AS 109][3].
///
/// [1]: http://people.sc.fsu.edu/~jburkardt/c_src/asa109/asa109.html
/// [2]: http://www.jstor.org/stable/2346798
/// [3]: http://www.jstor.org/stable/2346887
pub fn inv_inc_beta(mut alpha: f64, mut p: f64, mut q: f64, ln_beta: f64) -> f64 {
    // Algorithm AS 64
    // http://www.jstor.org/stable/2346798
    //
    // An approximation x₀ to x if found from (cf. Scheffé and Tukey, 1944)
    //
    // 1 + x₀   4p + 2q - 2
    // ------ = -----------
    // 1 - x₀      χ²(α)
    //
    // where χ²(α) is the upper α point of the χ² distribution with 2q degrees
    // of freedom and is obtained from Wilson and Hilferty’s approximation (cf.
    // Wilson and Hilferty, 1931)
    //
    // χ²(α) = 2q (1 - 1/(9q) + y(α) sqrt(1/(9q)))^3,
    //
    // y(α) being Hastings’ approximation (cf. Hastings, 1955) for the upper α
    // point of the standard normal distribution. If χ²(α) < 0, then
    //
    // x₀ = 1 - ((1 - α)q B(p, q))^(1/q).
    //
    // Again if (4p + 2q - 2)/χ²(α) does not exceed 1, x₀ is obtained from
    //
    // x₀ = (αp B(p, q))^(1/p).
    //
    // The final solution is obtained by the Newton–Raphson method from the
    // relation
    //
    //                  f(x[i-1])
    // x[i] = x[i-1] - ----------
    //                 f'(x[i-1])
    //
    // where
    //
    // f(x) = I(x, p, q) - α.

    use std::num::Float;

    #[inline(always)]
    fn exp(x: f64) -> f64 { x.exp() }
    #[inline(always)]
    fn ln(x: f64) -> f64 { x.ln() }
    #[inline(always)]
    fn pow(x: f64, y: f64) -> f64 { x.powf(y) }
    #[inline(always)]
    fn pow10(x: i32) -> f64 { 10f64.powi(x) }
    #[inline(always)]
    fn sqrt(x: f64) -> f64 { x.sqrt() }

    // Remark AS R83
    // http://www.jstor.org/stable/2347779
    const SAE: i32 = -30;
    const FPU: f64 = 1e-30; // 10^SAE

    if alpha <= 0.0 {
        return 0.0;
    }
    if 1.0 <= alpha {
        return 1.0;
    }

    let mut x;
    let mut y;

    let flip = 0.5 < alpha;
    if flip {
        x = p;
        p = q;
        q = x;
        alpha = 1.0 - alpha;
    }

    x = sqrt(-ln(alpha * alpha));
    y = x - (2.30753 + 0.27061 * x) / (1.0 + (0.99229 + 0.04481 * x) * x);

    if 1.0 < p && 1.0 < q {
        // Remark AS R19 and Algorithm AS 109
        // http://www.jstor.org/stable/2346887
        //
        // For p and q > 1, the approximation given by Carter (1947), which
        // improves the Fisher–Cochran formula, is generally better. For other
        // values of p and q en empirical investigation has shown that the
        // approximation given in AS 64 is adequate.
        let r = (y * y - 3.0) / 6.0;
        let s = 1.0 / (2.0 * p - 1.0);
        let t = 1.0 / (2.0 * q - 1.0);
        let h = 2.0 / (s + t);
        let w = y * sqrt(h + r) / h - (t - s) * (r + 5.0 / 6.0 - 2.0 / (3.0 * h));
        x = p / (p + q * exp(2.0 * w));
    } else {
        let mut t = 1.0 / (9.0 * q);
        t = 2.0 * q * pow(1.0 - t + y * sqrt(t), 3.0);
        if t <= 0.0 {
            x = 1.0 - exp((ln((1.0 - alpha) * q) + ln_beta) / q);
        } else {
            t = 2.0 * (2.0 * p + q - 1.0) / t;
            if t <= 1.0 {
                x = exp((ln(alpha * p) + ln_beta) / p);
            } else {
                x = 1.0 - 2.0 / (t + 1.0);
            }
        }
    }

    if x < 0.0001 {
        x = 0.0001;
    } else if 0.9999 < x {
        x = 0.9999;
    }

    // Remark AS R83
    // http://www.jstor.org/stable/2347779
    let e = (-5.0 / p / p - 1.0 / pow(alpha, 0.2) - 13.0) as i32;
    let acu = if e > SAE { pow10(e) } else { FPU };

    let mut tx;
    let mut yprev = 0.0;
    let mut sq = 1.0;
    let mut prev = 1.0;

    'outer: loop {
        // Remark AS R19 and Algorithm AS 109
        // http://www.jstor.org/stable/2346887
        y = inc_beta(x, p, q, ln_beta);
        y = (y - alpha) * exp(ln_beta + (1.0 - p) * ln(x) + (1.0 - q) * ln(1.0 - x));

        // Remark AS R83
        // http://www.jstor.org/stable/2347779
        if y * yprev <= 0.0 {
            prev = if sq > FPU { sq } else { FPU };
        }

        // Remark AS R19 and Algorithm AS 109
        // http://www.jstor.org/stable/2346887
        let mut g = 1.0;
        loop {
            loop {
                let adj = g * y;
                sq = adj * adj;

                if sq < prev {
                    tx = x - adj;
                    if 0.0 <= tx && tx <= 1.0 {
                        break;
                    }
                }
                g /= 3.0;
            }

            if prev <= acu || y * y <= acu {
                x = tx;
                break 'outer;
            }

            if tx != 0.0 && tx != 1.0 {
                break;
            }

            g /= 3.0;
        }

        if tx == x {
            break;
        }

        x = tx;
        yprev = y;
    }

    if flip { 1.0 - x } else { x }
}

#[cfg(test)]
mod test {
    #[phase(plugin)] extern crate assert;

    #[test]
    fn ln_beta() {
        let x = vec![0.25, 0.50, 0.75, 1.00];
        let y = vec![0.50, 0.75, 1.00, 1.25];
        let z = vec![
            1.6571065161914820,  0.8739177307778084,
            0.2876820724517809, -0.2231435513142098,
        ];

        assert_close!(x.iter().zip(y.iter()).map(|(&x, &y)| {
            super::ln_beta(x, y)
        }).collect::<Vec<_>>(), z);
    }

    #[test]
    fn inc_beta_small() {
        let (p, q) = (0.1, 0.2);
        let ln_beta = super::ln_beta(p, q);

        let x = vec![
            0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
            0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00,
        ];
        let y = vec![
            0.000000000000000e+00, 5.095391215346399e-01, 5.482400859052436e-01,
            5.732625733722232e-01, 5.925346573554778e-01, 6.086596697678208e-01,
            6.228433547203172e-01, 6.357578563479236e-01, 6.478288604374864e-01,
            6.593557133297501e-01, 6.705707961028990e-01, 6.816739425887479e-01,
            6.928567823206671e-01, 7.043251807250750e-01, 7.163269829958610e-01,
            7.291961263917867e-01, 7.434379555965913e-01, 7.599272566076309e-01,
            7.804880320024465e-01, 8.104335200313719e-01, 1.000000000000000e+00,
        ];

        assert_close!(x.iter().map(|&x| {
           super::inc_beta(x, p, q, ln_beta)
        }).collect::<Vec<_>>(), y);
    }

    #[test]
    fn inc_beta_large() {
        let (p, q) = (2.0, 3.0);
        let ln_beta = super::ln_beta(p, q);

        let x = vec![
            0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
            0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00,
        ];
        let y = vec![
            0.000000000000000e+00, 1.401875000000000e-02, 5.230000000000002e-02,
            1.095187500000000e-01, 1.807999999999999e-01, 2.617187500000001e-01,
            3.483000000000000e-01, 4.370187500000001e-01, 5.248000000000003e-01,
            6.090187500000001e-01, 6.875000000000000e-01, 7.585187500000001e-01,
            8.208000000000000e-01, 8.735187499999999e-01, 9.163000000000000e-01,
            9.492187500000000e-01, 9.728000000000000e-01, 9.880187500000001e-01,
            9.963000000000000e-01, 9.995187500000000e-01, 1.000000000000000e+00,
        ];

        assert_close!(x.iter().map(|&x| {
            super::inc_beta(x, p, q, ln_beta)
        }).collect::<Vec<_>>(), y);
    }

    #[test]
    fn inv_inc_beta_small() {
        let (p, q) = (0.2, 0.3);
        let ln_beta = super::ln_beta(p, q);

        let x = vec![
            0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
            0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00,
        ];
        let y = vec![
            0.000000000000000e+00, 2.793072851850660e-06, 8.937381711316164e-05,
            6.784491773826951e-04, 2.855345858289119e-03, 8.684107512129325e-03,
            2.144658503798324e-02, 4.568556852983932e-02, 8.683942933344659e-02,
            1.502095712585510e-01, 2.391350361479824e-01, 3.527066234122371e-01,
            4.840600731467657e-01, 6.206841200371190e-01, 7.474718280552188e-01,
            8.514539745840592e-01, 9.257428898178934e-01, 9.707021084050310e-01,
            9.923134416335146e-01, 9.992341305241808e-01, 1.000000000000000e+00,
        ];

        assert_close!(x.iter().map(|&x| {
            super::inv_inc_beta(x, p, q, ln_beta)
        }).collect::<Vec<_>>(), y);
    }

    #[test]
    fn inv_inc_beta_large() {
        let (p, q) = (1.0, 2.0);
        let ln_beta = super::ln_beta(p, q);

        let x = vec![
            0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
            0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00,
        ];
        let y = vec![
            0.000000000000000e+00, 0.025320565519104e+00, 0.051316701949486e+00,
            0.078045554270711e+00, 0.105572809000084e+00, 0.133974596215561e+00,
            0.163339973465924e+00, 0.193774225170145e+00, 0.225403330758517e+00,
            0.258380151290432e+00, 0.292893218813452e+00, 0.329179606750063e+00,
            0.367544467966324e+00, 0.408392021690038e+00, 0.452277442494834e+00,
            0.500000000000000e+00, 0.552786404500042e+00, 0.612701665379257e+00,
            0.683772233983162e+00, 0.776393202250021e+00, 1.000000000000000e+00,
        ];

        assert_close!(x.iter().map(|&x| {
            super::inv_inc_beta(x, p, q, ln_beta)
        }).collect::<Vec<_>>(), y);
    }
}

#[cfg(test)]
mod bench {
    extern crate test;

    use std::rand::random;

    #[bench]
    fn inc_beta(b: &mut test::Bencher) {
        let (p, q) = (0.5, 1.5);
        let ln_beta = super::ln_beta(p, q);
        let x = range(0u, 1000).map(|_| random()).collect::<Vec<_>>();

        b.iter(|| {
            for &x in x.iter() {
                test::black_box(super::inc_beta(x, p, q, ln_beta))
            }
        });
    }

    #[bench]
    fn inv_inc_beta(b: &mut test::Bencher) {
        let (p, q) = (0.5, 1.5);
        let ln_beta = super::ln_beta(p, q);
        let x = range(0u, 1000).map(|_| random()).collect::<Vec<_>>();

        b.iter(|| {
            for &x in x.iter() {
                test::black_box(super::inv_inc_beta(x, p, q, ln_beta))
            }
        });
    }
}
