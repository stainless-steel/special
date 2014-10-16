/// Computes the natural logarithm of the beta function.
pub fn log_beta(x: f64, y: f64) -> f64 {
    use super::{log_gamma};
    log_gamma(x) + log_gamma(y) - log_gamma(x + y)
}

/// Computes the incomplete beta function.
///
/// The code is based on a [C implementation][1] of the incomplete beta
/// function by [John Burkardt][2]. The original algorithm was published in
/// Applied Statistics and is known as [Algorithm AS 63][3]. The algorithm is
/// outlined in what follows.
///
/// The function uses the method discussed by Soper (1921). If p is not less
/// than (p + q)x and the integral part of q + (1 - x)(p + q) is a positive
/// integer, say s, reductions are made up to s times “by parts” using the
/// recurrence relation
///
/// ```ignore
///                 Γ(p+q)
/// I(x, p, q) = ----------- x^p (1-x)^(q-1) + I(x, p+1, q-1)
///              Γ(p+1) Γ(q)
/// ```
///
/// and then reductions are continued by “raising p” with the recurrence
/// relation
///
/// ```ignore
///                      Γ(p+q)
/// I(x, p+s, q-s) = --------------- x^(p+s) (1 - x)^(q-s) + I(x, p+s+1, q-s)
///                  Γ(p+s+1) Γ(q-s)
/// ```
///
/// If s is not a positive integer, reductions are made only by “raising p.”
/// The process of reduction is terminated when the relative contribution to the
/// integral is not greater than the value of ACU. If p is less than (p + q)x,
/// I(1-x, q, p) is first calculated by the above procedure and then I(x, p, q)
/// is obtained from the relation
///
/// ```ignore
/// I(x, p, q) = 1 - I(1-x, p, q).
/// ```
///
/// Soper (1921) demonstrated that the expansion of I(x, p, q) by “parts” and
/// “raising p” method as described above converges more rapidly than any other
/// series expansions.
///
/// [1]: http://people.sc.fsu.edu/~jburkardt/c_src/asa109/asa109.html
/// [2]: http://people.sc.fsu.edu/~jburkardt/i.html
/// [3]: http://www.jstor.org/stable/2346797
pub fn inc_beta(x: f64, mut p: f64, mut q: f64, log_beta: f64) -> f64 {
    const ACU: f64 = 0.1e-14;

    if x <= 0.0 { return 0.0; }
    if 1.0 <= x { return 1.0; }

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

        temp = ::std::num::abs(term);
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
    alpha = unsafe {
        use super::m::{log, exp};
        alpha * exp(p * log(pbase) + (q - 1.0) * log(qbase) - log_beta) / p
    };

    if flip { 1.0 - alpha } else { alpha }
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

    #[test]
    fn log_beta() {
        let x = vec![0.25, 0.50, 0.75, 1.00];
        let y = vec![0.50, 0.75, 1.00, 1.25];
        let z = vec![1.6571065161914820, 0.8739177307778084,
                     0.2876820724517809, -0.2231435513142098];

        assert_close!(x.iter().zip(y.iter()).map(|(&x, &y)| {
            super::log_beta(x, y)
        }).collect(), z);
    }

    #[test]
    fn inc_beta_small() {
        let (p, q) = (0.1, 0.2);
        let log_beta = super::log_beta(p, q);

        let x = vec![
            0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45,
            0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00,
        ];
        let y = vec![
            0.000000000000000e+00, 5.095391215346399e-01,
            5.482400859052436e-01, 5.732625733722232e-01,
            5.925346573554778e-01, 6.086596697678208e-01,
            6.228433547203172e-01, 6.357578563479236e-01,
            6.478288604374864e-01, 6.593557133297501e-01,
            6.816739425887479e-01, 6.928567823206671e-01,
            7.043251807250750e-01, 7.163269829958610e-01,
            7.291961263917867e-01, 7.434379555965913e-01,
            7.599272566076309e-01, 7.804880320024465e-01,
            8.104335200313719e-01, 1.000000000000000e+00,
        ];

        assert_close!(x.iter().map(|&x| {
           super::inc_beta(x, p, q, log_beta)
        }).collect(), y);
    }

    #[test]
    fn inc_beta_large() {
        let (p, q) = (2.0, 3.0);
        let log_beta = super::log_beta(p, q);

        let x = vec![
            0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45,
            0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00,
        ];
        let y = vec![
            0.000000000000000e+00, 1.401875000000000e-02,
            5.230000000000002e-02, 1.095187500000000e-01,
            1.807999999999999e-01, 2.617187500000001e-01,
            3.483000000000000e-01, 4.370187500000001e-01,
            5.248000000000003e-01, 6.090187500000001e-01,
            7.585187500000001e-01, 8.208000000000000e-01,
            8.735187499999999e-01, 9.163000000000000e-01,
            9.492187500000000e-01, 9.728000000000000e-01,
            9.880187500000001e-01, 9.963000000000000e-01,
            9.995187500000000e-01, 1.000000000000000e+00,
        ];

        assert_close!(x.iter().map(|&x| {
            super::inc_beta(x, p, q, log_beta)
        }).collect(), y);
    }
}
