/// Lambert W function.
pub trait Omega
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
    /// use special::Omega;
    ///
    /// assert!((1.0.lambert_w0() - 0.5671432904097838).abs() < 1e-15);
    /// ```
    ///
    /// ## References
    ///
    /// - T. Fukushima, Precise and fast computation of Lambert W function
    /// by piecewise minimax rational function approximation with variable transformation,
    /// To be submitted, 2020, <https://doi.org/10.13140/RG.2.2.30264.37128>.  
    /// - The [`lambert_w`] crate.
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
    /// use special::Omega;
    ///
    /// assert!(((-f64::ln(2.0) / 2.0).lambert_wm1() + f64::ln(4.0)).abs() < 1e-15);
    /// ```
    ///
    /// ## References
    ///
    /// - T. Fukushima, Precise and fast computation of Lambert W function
    /// by piecewise minimax rational function approximation with variable transformation,
    /// To be submitted, 2020, <https://doi.org/10.13140/RG.2.2.30264.37128>.  
    /// - The [`lambert_w`] crate.
    fn lambert_wm1(self) -> Self;
}

impl Omega for f64 {
    fn lambert_w0(self) -> Self {
        lambert_w::lambert_w0(self)
    }

    fn lambert_wm1(self) -> Self {
        lambert_w::lambert_wm1(self)
    }
}
