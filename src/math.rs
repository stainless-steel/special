macro_rules! declare {
    ($(($method:ident, $f32:expr, $f64:expr, ($($argument:ident: $type_outer:tt as $type_inner:tt),*) -> $return:tt),)*) => {
        pub trait Float: Sized {
            const PI: Self;

            $(
                fn $method(self, $($argument: $type_outer),*) -> $return;
            )*
        }
    }
}

macro_rules! implement {
    ($(($method:ident, $f32:expr, $f64:expr, ($($argument:ident: $type_outer:tt as $type_inner:tt),*) -> $return:tt),)*) => {
        impl Float for f32 {
            const PI: Self = core::f32::consts::PI;

            $(
                #[inline(always)]
                fn $method(self, $($argument: $type_outer),*) -> $return {
                    $f32(self, $($argument as $type_inner),*)
                }
            )*
        }

        impl Float for f64 {
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
            (erf, libm::erff, libm::erf, () -> Self),
            (erfc, libm::erfcf, libm::erfc, () -> Self),
            (exp, libm::expf, libm::exp, () -> Self),
            (floor, libm::floorf, libm::floor, () -> Self),
            (lgamma, libm::lgammaf_r, libm::lgamma_r, () -> (Self, i32)),
            (ln, libm::logf, libm::log, () -> Self),
            (powf, libm::powf, libm::pow, (exponent: Self as Self) -> Self),
            (powi, libm::powf, libm::pow, (exponent: i32 as Self) -> Self),
            (sin, libm::sinf, libm::sin, () -> Self),
            (sqrt, libm::sqrtf, libm::sqrt, () -> Self),
            (tgamma, libm::tgammaf, libm::tgamma, () -> Self),
        }
    };
}

run! { declare }
run! { implement }
