macro_rules! declare {
    ($(($method:ident, $f32:expr, $f64:expr, ($($argument:ident: $type_outer:tt as $type_inner:tt),*) -> $return:tt),)*) => {
        /// Primitive functions.
        pub trait Primitive: Sized {
            /// The number π.
            const PI: Self;

            $(
                /// See `libm` for further details.
                fn $method(self, $($argument: $type_outer),*) -> $return;
            )*
        }
    }
}

macro_rules! implement {
    ($(($method:ident, $f32:expr, $f64:expr, ($($argument:ident: $type_outer:tt as $type_inner:tt),*) -> $return:tt),)*) => {
        impl Primitive for f32 {
            const PI: Self = core::f32::consts::PI;

            $(
                #[inline(always)]
                fn $method(self, $($argument: $type_outer),*) -> $return {
                    $f32(self, $($argument as $type_inner),*)
                }
            )*
        }

        impl Primitive for f64 {
            const PI: Self = core::f64::consts::PI;

            $(
                #[inline(always)]
                fn $method(self, $($argument: $type_outer),*) -> $return {
                    $f64(self, $($argument as $type_inner),*)
                }
            )*
        }
    };
}

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
