/// Lambert W functions.
pub trait LambertW
where
    Self: Sized,
{
    /// Compute the real-valued parts of the princial branch of the Lambert W function.
    ///
    /// This function is the inverse of the function
    ///
    /// ```math
    /// f(w) = we^w
    /// ```
    ///
    /// when w >= -1.
    ///
    /// ## Example
    ///
    /// ```
    /// use special::LambertW;
    ///
    /// let Ω = 0.5671432904097838;
    /// assert!((1.0.lambert_w0() - Ω).abs() < 1e-15);
    /// ```
    ///
    /// This code is based on a [Rust implementation][1] by Johanna Sörngård.
    ///
    /// ## References
    ///
    /// 1. T. Fukushima, Precise and fast computation of Lambert W function by piecewise minimax
    ///    rational function approximation with variable transformation.
    ///
    /// [1]: https://crates.io/crates/lambert_w
    fn lambert_w0(self) -> Self;

    /// Compute the real-valued parts of the secondary branch of the Lambert W function.
    ///
    /// This function is the inverse of the function
    ///
    /// ```math
    /// f(w) = we^w
    /// ```
    ///
    /// when w < -1.
    ///
    /// ## Example
    ///
    /// ```
    /// use special::LambertW;
    ///
    /// assert!(((-f64::ln(2.0) / 2.0).lambert_wm1() + f64::ln(4.0)).abs() < 1e-15);
    /// ```
    ///
    /// This code is based on a [Rust implementation][1] by Johanna Sörngård.
    ///
    /// ## References
    ///
    /// 1. T. Fukushima, Precise and fast computation of Lambert W function by piecewise minimax
    ///    rational function approximation with variable transformation.
    ///
    /// [1]: https://crates.io/crates/lambert_w
    fn lambert_wm1(self) -> Self;
}

impl LambertW for f64 {
    #[inline]
    fn lambert_w0(self) -> Self {
        lambert_w::lambert_w0(self)
    }

    #[inline]
    fn lambert_wm1(self) -> Self {
        lambert_w::lambert_wm1(self)
    }
}
