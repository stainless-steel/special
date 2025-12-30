macro_rules! declare_method {
    ($(#[$attribute:meta])* $name:ident($($argument:ident),*)) => {
        $(#[$attribute])*
        fn $name(self, $($argument: Self,)*) -> Self;
    };
}

macro_rules! define_method {
    (@first_m $name:ident -> $backend:ident($($argument:ident),*)) => {
        #[inline]
        fn $name(self, $($argument: Self,)*) -> Self {
            ellip::$backend(self, $($argument,)*).unwrap()
        }
    };
    (@first_kc $name:ident -> $backend:ident($($argument:ident),*)) => {
        #[inline]
        fn $name(self, $($argument: Self,)*) -> Self {
            debug_assert!(self > 0.0, "m (self) cannot be less than 1");
            ellip::$backend((1.0 - self).sqrt(), $($argument,)*).unwrap()
        }
    };
    (@second_kc $name:ident -> $backend:ident($x:ident $(, $argument:ident)*)) => {
        #[inline]
        fn $name(self, $x: Self, $($argument: Self,)*) -> Self {
            debug_assert!(self > 0.0, "m (self) cannot be less than 1");
            ellip::$backend($x, (1.0 - self).sqrt(), $($argument,)*).unwrap()
        }
    };
    ($name:ident -> $backend:ident($($argument:ident),*)) => {
        #[inline]
        fn $name(self, $($argument: Self,)*) -> Self {
            ellip::$backend($($argument,)* self).unwrap()
        }
    };
}

macro_rules! implement {
    ($(
        $(#[$attribute:meta])*
        $(@$flag:ident)? $name:ident -> $backend:ident($($argument:tt)*),
    )*) => {
        /// Elliptic integrals.
        pub trait Elliptic: Sized {
            $(declare_method!($(#[$attribute])* $name($($argument)*));)*
        }

        impl Elliptic for f32 {
            $(define_method!($(@$flag)? $name -> $backend($($argument)*));)*
        }

        impl Elliptic for f64 {
            $(define_method!($(@$flag)? $name -> $backend($($argument)*));)*
        }
    };
}

implement!(
    /// Compute the complete elliptic integral in Legendre's form of the first kind (K).
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert.
    ///
    /// ## Parameters
    ///
    /// - `self`: elliptic parameter (m)
    ///
    /// ## Examples
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let m = 0.5;
    /// assert::close(m.elliptic_k(), 1.8540746773013719, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if m > 1.
    ///
    /// [1]: https://crates.io/crates/ellip
    elliptic_k -> ellipk(),

    /// Compute the complete elliptic integral in Legendre's form of the second kind (E).
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert.
    ///
    /// ## Parameters
    ///
    /// - `self`: elliptic parameter (m)
    ///
    /// ## Examples
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let m = 0.5;
    /// assert::close(m.elliptic_e(), 1.3506438810476755, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if m > 1.
    ///
    /// [1]: https://crates.io/crates/ellip
    elliptic_e -> ellipe(),

    /// Compute the complete elliptic integral in Legendre's form of the third kind (Π).
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert.
    ///
    /// ## Parameters
    ///
    /// - `self`: elliptic parameter (m)
    /// - `n`: characteristic (n)
    ///
    /// ## Examples
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let m = 0.5;
    /// let n = 0.5;
    /// assert::close(m.elliptic_pi(n), 2.7012877620953506, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if m > 1.
    ///
    /// [1]: https://crates.io/crates/ellip
    elliptic_pi -> ellippi(n),

    /// Compute the complete elliptic integral in Legendre's form of Legendre's type (D).
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert.
    ///
    /// ## Parameters
    ///
    /// - `self`: elliptic parameter (m)
    ///
    /// ## Examples
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let m = 0.5;
    /// assert::close(m.elliptic_d(), 1.0068615925073927, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if m > 1.
    ///
    /// [1]: https://crates.io/crates/ellip
    elliptic_d -> ellipd(),

    /// Compute the incomplete elliptic integral in Legendre's form of the first kind (F).
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert.
    ///
    /// ## Parameters
    ///
    /// - `self`: elliptic parameter (m)
    /// - `phi`: amplitude angle (φ)
    ///
    /// ## Examples
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let m = 0.5;
    /// let phi = std::f64::consts::FRAC_PI_4;
    /// assert::close(m.inc_elliptic_f(phi), 0.826017876249245, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if m sin²φ > 1.
    ///
    /// [1]: https://crates.io/crates/ellip
    inc_elliptic_f -> ellipf(phi),

    /// Compute the incomplete elliptic integral in Legendre's form of the second kind (E).
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert.
    ///
    /// ## Parameters
    ///
    /// - `self`: elliptic parameter (m)
    /// - `phi`: amplitude angle (φ)
    ///
    /// ## Examples
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let m = 0.5;
    /// let phi = std::f64::consts::FRAC_PI_4;
    /// assert::close(m.inc_elliptic_e(phi), 0.7481865041776612, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if m sin²φ > 1.
    ///
    /// [1]: https://crates.io/crates/ellip
    inc_elliptic_e -> ellipeinc(phi),

    /// Compute the incomplete elliptic integral in Legendre's form of the third kind (Π).
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert.
    ///
    /// ## Parameters
    ///
    /// - `self`: elliptic parameter (m)
    /// - `phi`: amplitude angle (φ)
    /// - `n`: characteristic
    ///
    /// ## Examples
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let m = 0.5;
    /// let phi = std::f64::consts::FRAC_PI_4;
    /// let n = 0.5;
    /// assert::close(m.inc_elliptic_pi(phi, n), 0.9190227391656969, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if m sin²φ > 1, n sin²φ = 1, or m ≥ 1 and φ is not a multiple of π/2.
    ///
    /// [1]: https://crates.io/crates/ellip
    inc_elliptic_pi -> ellippiinc(phi, n),

    /// Compute the incomplete elliptic integral in Legendre's form of the third kind (Π) using
    /// Bulirsch's method.
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert.
    ///
    /// ## Parameters
    ///
    /// - `self`: elliptic parameter (m)
    /// - `phi`: amplitude angle (φ)
    /// - `n`: characteristic
    ///
    /// ## Examples
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let m = 0.5;
    /// let phi = std::f64::consts::FRAC_PI_4;
    /// let n = 0.5;
    /// assert::close(m.inc_elliptic_pi_bulirsch(phi, n), 0.9190227391656969, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if m sin²φ > 1, n sin²φ = 1, or m ≥ 1 and φ is not a multiple of π/2.
    ///
    /// [1]: https://crates.io/crates/ellip
    inc_elliptic_pi_bulirsch -> ellippiinc_bulirsch(phi, n),

    /// Compute the incomplete elliptic integral in Legendre's form of Legendre's type (D).
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert.
    ///
    /// ## Parameters
    ///
    /// - `self`: elliptic parameter (m)
    /// - `phi`: amplitude angle (φ)
    ///
    /// ## Examples
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let m = 0.5;
    /// let phi = std::f64::consts::FRAC_PI_4;
    /// assert::close(m.inc_elliptic_d(phi), 0.15566274414316758, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if m sin²φ > 1.
    ///
    /// [1]: https://crates.io/crates/ellip
    inc_elliptic_d -> ellipdinc(phi),

    /// Compute the general complete elliptic integral in Bulirsch's form.
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert. Note that the original
    /// literature by Bulirsch used the complementary modulus kc where kc = √(1 - m).
    ///
    /// ## Parameters
    ///
    /// - `self`: elliptic parameter (m)
    /// - `p`: parameter
    /// - `a` and `b`: coefficients
    ///
    /// ## Examples
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let m = 0.5;
    /// let p = 1.0;
    /// let a = 1.0;
    /// let b = 1.0;
    ///
    /// assert::close(m.elliptic_bulirsch(p, a, b), 1.8540746773013717, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if m = 1, p = 0, or more than one arguments are infinite.
    ///
    /// [1]: https://crates.io/crates/ellip
    @first_kc elliptic_bulirsch -> cel(p, a, b),

    /// Compute the complete elliptic integral in Bulirsch's form of the first kind.
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert. Note that the original
    /// literature by Bulirsch used the complementary modulus kc where kc = √(1 - m).
    ///
    /// ## Parameters
    ///
    /// - `self`: elliptic parameter (m)
    ///
    /// ## Examples
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let m = 0.5;
    ///
    /// assert::close(m.elliptic_bulirsch_1(), 1.8540746773013717, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if m = 1.
    ///
    /// [1]: https://crates.io/crates/ellip
    @first_kc elliptic_bulirsch_1 -> cel1(),

    /// Compute the complete elliptic integral in Bulirsch's form of the second kind.
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert. Note that the original
    /// literature by Bulirsch used the complementary modulus kc where kc = √(1 - m).
    ///
    /// ## Parameters
    ///
    /// - `self`: elliptic parameter (m)
    /// - `a` and `b`: coefficients
    ///
    /// ## Examples
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let m = 0.5;
    /// let a = 1.0;
    /// let b = 1.0;
    ///
    /// assert::close(m.elliptic_bulirsch_2(a, b), 1.8540746773013717, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if m = 1 or more than one arguments are infinite.
    ///
    /// [1]: https://crates.io/crates/ellip
    @first_kc elliptic_bulirsch_2 -> cel2(a, b),

    /// Compute the incomplete elliptic integral Bulirsch's form of the first kind.
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert. Note that the original
    /// literature by Bulirsch used the complementary modulus kc where kc = √(1 - m).
    ///
    /// ## Parameters
    ///
    /// - `self`: elliptic parameter (m)
    /// - `x`: tangent of the amplitude angle
    ///
    /// ## Examples
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let m = 0.5;
    /// let x = std::f64::consts::FRAC_PI_4.tan();
    ///
    /// assert::close(m.inc_elliptic_bulirsch_1(x), 0.826017876249245, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if m = 1.
    ///
    /// [1]: https://crates.io/crates/ellip
    @second_kc inc_elliptic_bulirsch_1 -> el1(x),

    /// Compute the incomplete elliptic integral in Bulirsch's form of the second kind.
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert. Note that the original
    /// literature by Bulirsch used the complementary modulus kc where kc = √(1 - m).
    ///
    /// ## Parameters
    ///
    /// - `self`: elliptic parameter (m)
    /// - `x`: tangent of the amplitude angle
    /// - `a` and `b`: coefficients
    ///
    /// ## Examples
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let m = 0.5;
    /// let x = std::f64::consts::FRAC_PI_4.tan();
    /// let a = 1.0;
    /// let b = 1.0;
    ///
    /// assert::close(m.inc_elliptic_bulirsch_2(x, a, b), 0.826017876249245, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if m = 1.
    ///
    /// [1]: https://crates.io/crates/ellip
    @second_kc inc_elliptic_bulirsch_2 -> el2(x, a, b),

    /// Compute the incomplete elliptic integral Bulirsch's form of the third kind.
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert. Note that the original
    /// literature by Bulirsch used the complementary modulus kc where kc = √(1 - m).
    ///
    /// ## Parameters
    ///
    /// - `self`: elliptic parameter (m)
    /// - `x`: tangent of the amplitude angle
    /// - `p`: a parameter
    ///
    /// ## Examples
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let m = 0.5;
    /// let x = std::f64::consts::FRAC_PI_4.tan();
    /// let p = 1.0;
    ///
    /// assert::close(m.inc_elliptic_bulirsch_3(x, p), 0.826017876249245, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if m = 1, 1 + px² = 0, or m < 0 for p < 0.
    ///
    /// [1]: https://crates.io/crates/ellip
    @second_kc inc_elliptic_bulirsch_3 -> el3(x, p),

    /// Compute the symmetric elliptic integral in Carlson's form of the first kind (RF).
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert.
    ///
    /// ## Parameters
    ///
    /// - `self`, `y`, and `z`: symmetric arguments
    ///
    /// ## Examples
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let x = 1.0;
    /// let y = 0.5;
    /// let z = 0.25;
    ///
    /// assert::close(x.elliptic_carlson_rf(y, z), 1.370171633266872, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if any of x, y, or z is negative, or more than one of them are zero.
    ///
    /// [1]: https://crates.io/crates/ellip
    elliptic_carlson_rf -> elliprf(y, z),

    /// Compute the symmetric elliptic integral in Carlson's form of the second kind (RG).
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert.
    ///
    /// ## Parameters
    ///
    /// - `self`, `y`, and `z`: symmetric arguments
    ///
    /// ## Examples
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let x = 1.0;
    /// let y = 0.5;
    /// let z = 0.25;
    ///
    /// assert::close(x.elliptic_carlson_rg(y, z), 0.7526721491833781, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if any of x, y, or z is negative or infinite.
    ///
    /// [1]: https://crates.io/crates/ellip
    elliptic_carlson_rg -> elliprg(y, z),

    /// Compute the symmetric elliptic integral in Carlson's form of the third kind (RJ).
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert.
    ///
    /// ## Parameters
    ///
    /// - `self`, `y`, and `z`: symmetric arguments
    /// - `p`: parameter
    ///
    /// ## Examples
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let x = 1.0;
    /// let y = 0.5;
    /// let z = 0.25;
    /// let p = 0.125;
    ///
    /// assert::close(x.elliptic_carlson_rj(y, z, p), 5.680557292035963, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if p = 0, any of x, y, or z is negative, or more than one of them are zero.
    ///
    /// [1]: https://crates.io/crates/ellip
    @first_m elliptic_carlson_rj -> elliprj(y, z, p),

    /// Compute the degenerate symmetric elliptic integral in Carlson's form (RC).
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert.
    ///
    /// ## Parameters
    ///
    /// - `self` and `y`: symmetric arguments
    ///
    /// ## Examples
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let x = 1.0;
    /// let y = 0.5;
    ///
    /// assert::close(x.elliptic_carlson_rc(y), 1.2464504802804608, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if x < 0, y = 0, or y < 0.
    ///
    /// [1]: https://crates.io/crates/ellip
    @first_m elliptic_carlson_rc -> elliprc(y),

    /// Compute the symmetric elliptic integral in Carlson's form of the second kind (RD).
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert.
    ///
    /// ## Parameters
    ///
    /// - `self`, `y`, and `z`: symmetric arguments
    ///
    /// ## Examples
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let x = 1.0;
    /// let y = 0.5;
    /// let z = 0.25;
    ///
    /// assert::close(x.elliptic_carlson_rd(y, z), 4.022594757168912, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if x < 0, y < 0, z ≤ 0, or when both x and y are zero.
    ///
    /// [1]: https://crates.io/crates/ellip
    @first_m elliptic_carlson_rd -> elliprd(y, z),

    /// Compute Jacobi's Zeta function.
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert.
    ///
    /// ## Parameters
    ///
    /// - `self`: elliptic parameter (m)
    /// - `phi`: amplitude angle (φ)
    ///
    /// ## Examples
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let m = 0.5;
    /// let phi = std::f64::consts::FRAC_PI_4;
    ///
    /// assert::close(m.jacobi_zeta(phi), 0.146454543836188, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if m > 1, φ is infinite, or m is infinite.
    ///
    /// [1]: https://crates.io/crates/ellip
    jacobi_zeta -> jacobi_zeta(phi),

    /// Compute Heuman's Lambda function.
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert.
    ///
    /// ## Parameters
    ///
    /// - `self`: elliptic parameter (m)
    /// - `phi`: amplitude angle (φ)
    ///
    /// ## Examples
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let m = 0.5;
    /// let phi = std::f64::consts::FRAC_PI_4;
    ///
    /// assert::close(m.heuman_lambda(phi), 0.6183811341833665, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if m < 0, m ≥ 1, or phi is infinite.
    ///
    /// [1]: https://crates.io/crates/ellip
    heuman_lambda -> heuman_lambda(phi),
);
