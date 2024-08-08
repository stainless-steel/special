macro_rules! declare {
    ($(($name:ident, $f32:expr, $f64:expr, ($($argument:ident: $type_outer:ty as $type_inner:ty),*) -> $return:ty),)*) => {
        /// Primitive functions.
        pub trait Primitive: Sized {
            const PI: Self;
            $(fn $name(self, $($argument: $type_outer),*) -> $return;)*
        }
    }
}

macro_rules! implement {
    ($(($name:ident, $f32:expr, $f64:expr, ($($argument:ident: $type_outer:ty as $type_inner:ty),*) -> $return:ty),)*) => {
        impl Primitive for f32 {
            const PI: Self = core::f32::consts::PI;
            $(implement! { @method $name, $f32, ($($argument: $type_outer as $type_inner),*) -> $return })*
        }

        impl Primitive for f64 {
            const PI: Self = core::f64::consts::PI;
            $(implement! { @method $name, $f64, ($($argument: $type_outer as $type_inner),*) -> $return })*
        }
    };
    (@method $name:ident, $function:expr, ($($argument:ident: $type_outer:ty as $type_inner:ty),*) -> $return:ty) => {
        #[inline(always)]
        fn $name(self, $($argument: $type_outer),*) -> $return {
            $function(self, $($argument as $type_inner),*)
        }
    }
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
