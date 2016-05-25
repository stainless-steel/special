use math;

/// Error functions.
pub trait Error {
    /// Compute the error function.
    fn erf(self) -> Self;

    /// Compute the complementary error function.
    fn erfc(self) -> Self;
}

macro_rules! implement {
    ($($kind:ty),*) => ($(
        impl Error for $kind {
            #[inline]
            fn erf(self) -> Self {
                unsafe { math::erf(self as f64) as Self }
            }

            #[inline]
            fn erfc(self) -> Self {
                unsafe { math::erfc(self as f64) as Self }
            }
        }
    )*);
}

implement!(f32, f64);
