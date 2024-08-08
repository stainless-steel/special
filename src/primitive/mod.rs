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

#[cfg_attr(feature = "no_std", path = "extrinsic.rs")]
#[cfg_attr(not(feature = "no_std"), path = "intrinsic.rs")]
mod implementation;

pub use implementation::Primitive;
