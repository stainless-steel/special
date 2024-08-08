macro_rules! run {
    ($macro: ident) => {
        $macro! {
            (abs, libm::fabsf, libm::fabs, () -> Self),
            (atan, libm::atanf, libm::atan, () -> Self),
            (erf, libm::erff, libm::erf, () -> Self),
            (erfc, libm::erfcf, libm::erfc, () -> Self),
            (exp, libm::expf, libm::exp, () -> Self),
            (exp_m1, libm::expm1f, libm::expm1, () -> Self),
            (floor, libm::floorf, libm::floor, () -> Self),
            (lgamma, libm::lgammaf_r, libm::lgamma_r, () -> (Self, i32)),
            (ln, libm::logf, libm::log, () -> Self),
            (ln_1p, libm::log1pf, libm::log1p, () -> Self),
            (powf, libm::powf, libm::pow, (exponent: Self as Self) -> Self),
            (powi, libm::powf, libm::pow, (exponent: i32 as Self) -> Self),
            (round, libm::roundf, libm::round, () -> Self),
            (sin, libm::sinf, libm::sin, () -> Self),
            (sqrt, libm::sqrtf, libm::sqrt, () -> Self),
            (tan, libm::tanf, libm::tan, () -> Self),
            (tgamma, libm::tgammaf, libm::tgamma, () -> Self),
            (trunc, libm::truncf, libm::trunc, () -> Self),
        }
    };
}

run! { declare }
run! { implement }
