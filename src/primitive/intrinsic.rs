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
