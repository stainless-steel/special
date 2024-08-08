macro_rules! run {
    ($macro: ident) => {
        $macro! {
            (abs, Self::abs, Self::abs, () -> Self),
            (atan, Self::atan, Self::atan, () -> Self),
            (erf, libm::erff, libm::erf, () -> Self),
            (erfc, libm::erfcf, libm::erfc, () -> Self),
            (exp, Self::exp, Self::exp, () -> Self),
            (exp_m1, Self::exp_m1, Self::exp_m1, () -> Self),
            (floor, Self::floor, Self::floor, () -> Self),
            (lgamma, libm::lgammaf_r, libm::lgamma_r, () -> (Self, i32)),
            (ln, Self::ln, Self::ln, () -> Self),
            (ln_1p, Self::ln_1p, Self::ln_1p, () -> Self),
            (powf, Self::powf, Self::powf, (exponent: Self as Self) -> Self),
            (powi, Self::powi, Self::powi, (exponent: i32 as i32) -> Self),
            (round, Self::round, Self::round, () -> Self),
            (sin, Self::sin, Self::sin, () -> Self),
            (sqrt, Self::sqrt, Self::sqrt, () -> Self),
            (tan, Self::tan, Self::tan, () -> Self),
            (tgamma, libm::tgammaf, libm::tgamma, () -> Self),
            (trunc, Self::trunc, Self::trunc, () -> Self),
        }
    };
}

run! { declare }
run! { implement }
