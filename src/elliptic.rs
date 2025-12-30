macro_rules! declare_method {
    ($(#[$attribute:meta])* $name:ident($($argument:ident),*)) => {
        $(#[$attribute])*
        fn $name(self, $($argument: Self,)*) -> Self;
    };
}

// By default, define_method expect self as the last argument.
// With @kc, define_method expect self as **m** as the first argument.
// With @kc_x, define_method expect x as the first and self as the second argument.
// With @first, define_method expect self as the first argument.
macro_rules! define_method {
    (@kc $name:ident -> $backend:ident($($argument:ident),*)) => {
        fn $name(self, $($argument: Self,)*) -> Self {
            debug_assert!(self > 0.0, concat!("m (self) cannot be less than 1"));
            ellip::$name((1.0 - self).sqrt(), $($argument,)*).unwrap()
        }
    };
    (@kc_x $name:ident -> $backend:ident($x:ident)) => {
        fn $name(self, $x: Self) -> Self {
            debug_assert!(self > 0.0, concat!("m (self) cannot be less than 1"));
            ellip::$name($x, (1.0 - self).sqrt()).unwrap()
        }
    };
    (@kc_x $name:ident -> $backend:ident($x:ident, $($argument:ident),+)) => {
        fn $name(self, $x: Self, $($argument: Self,)*) -> Self {
            debug_assert!(self > 0.0, concat!("m (self) cannot be less than 1"));
            ellip::$name($x, (1.0 - self).sqrt(), $($argument,)*).unwrap()
        }
    };
    (@first $name:ident -> $backend:ident($($argument:ident),*)) => {
        fn $name(self, $($argument: Self,)*) -> Self {
            ellip::$name(self, $($argument,)*).unwrap()
        }
    };
    ($name:ident -> $backend:ident($($argument:ident),*)) => {
        fn $name(self, $($argument: Self,)*) -> Self {
            ellip::$name($($argument,)* self).unwrap()
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
    /// Compute the complete elliptic integral of the first kind (K).
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert.
    ///
    /// ## Parameters
    ///
    /// - `self`: elliptic parameter (m)
    ///
    /// ## Example
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let m = 0.5;
    /// assert::close(m.ellipk(), 1.8540746773013719, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if m > 1.
    ///
    /// [1]: https://crates.io/crates/ellip
    ellipk -> ellipk(),

    /// Compute the complete elliptic integral of the second kind (E).
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert.
    ///
    /// ## Parameters
    ///
    /// - `self`: elliptic parameter (m)
    ///
    /// ## Example
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let m = 0.5;
    /// assert::close(m.ellipe(), 1.3506438810476755, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if m > 1.
    ///
    /// [1]: https://crates.io/crates/ellip
    ellipe -> ellipe(),

    /// Compute the complete elliptic integral of the third kind (Π).
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert.
    ///
    /// ## Parameters
    ///
    /// - `self`: elliptic parameter (m)
    /// - n: characteristic (n)
    ///
    /// ## Example
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let m = 0.5;
    /// let n = 0.5;
    /// assert::close(m.ellippi(n), 2.7012877620953506, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if m > 1.
    ///
    /// [1]: https://crates.io/crates/ellip
    ellippi -> ellippi(n),

    /// Compute the complete elliptic integral of Legendre's type (D).
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert.
    ///
    /// ## Parameters
    ///
    /// - `self`: elliptic parameter (m)
    ///
    /// ## Example
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let m = 0.5;
    /// assert::close(m.ellipd(), 1.0068615925073927, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if m > 1.
    ///
    /// [1]: https://crates.io/crates/ellip
    ellipd -> ellipd(),

    // <--- Incomplete Legendre's Integrals --->

    /// Compute the incomplete elliptic integral of the first kind (F).
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert.
    ///
    /// ## Parameters
    ///
    /// - `self`: elliptic parameter (m)
    /// - `phi`: amplitude angle (φ)
    ///
    /// ## Example
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let m = 0.5;
    /// let phi = std::f64::consts::FRAC_PI_4;
    /// assert::close(m.ellipf(phi), 0.826017876249245, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if m sin²φ > 1.
    ///
    /// [1]: https://crates.io/crates/ellip
    ellipf -> ellipf(phi),

    /// Compute the incomplete elliptic integral of the second kind (E).
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert.
    ///
    /// ## Parameters
    ///
    /// - `self`: elliptic parameter (m)
    /// - `phi`: amplitude angle (φ)
    ///
    /// ## Example
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let m = 0.5;
    /// let phi = std::f64::consts::FRAC_PI_4;
    /// assert::close(m.ellipeinc(phi), 0.7481865041776612, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if m sin²φ > 1.
    ///
    /// [1]: https://crates.io/crates/ellip
    ellipeinc -> ellipeinc(phi),

    /// Compute the incomplete elliptic integral of the third kind (Π).
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert.
    ///
    /// ## Parameters
    ///
    /// - `self`: elliptic parameter (m)
    /// - `phi`: amplitude angle (φ)
    /// - `n`: characteristic
    ///
    /// ## Example
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let m = 0.5;
    /// let phi = std::f64::consts::FRAC_PI_4;
    /// let n = 0.5;
    /// assert::close(m.ellippiinc(phi, n), 0.9190227391656969, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if m sin²φ > 1, n sin²φ = 1, or m ≥ 1 and φ is not a multiple of π/2.
    ///
    /// [1]: https://crates.io/crates/ellip
    ellippiinc -> ellippiinc(phi, n),

    /// Compute the incomplete elliptic integral of the third kind (Π) using Bulirsch's method.
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert.
    ///
    /// ## Parameters
    ///
    /// - `self`: elliptic parameter (m)
    /// - `phi`: amplitude angle (φ)
    /// - `n`: characteristic
    ///
    /// ## Example
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let m = 0.5;
    /// let phi = std::f64::consts::FRAC_PI_4;
    /// let n = 0.5;
    /// assert::close(m.ellippiinc_bulirsch(phi, n), 0.9190227391656969, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if m sin²φ > 1, n sin²φ = 1, or m ≥ 1 and φ is not a multiple of π/2.
    ///
    /// [1]: https://crates.io/crates/ellip
    ellippiinc_bulirsch -> ellippiinc_bulirsch(phi, n),

    /// Compute the incomplete elliptic integral of Legendre's type (D).
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert.
    ///
    /// ## Parameters
    ///
    /// - `self`: elliptic parameter (m)
    /// - `phi`: amplitude angle (φ)
    ///
    /// ## Example
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let m = 0.5;
    /// let phi = std::f64::consts::FRAC_PI_4;
    /// assert::close(m.ellipdinc(phi), 0.15566274414316758, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if m sin²φ > 1.
    ///
    /// [1]: https://crates.io/crates/ellip
    ellipdinc -> ellipdinc(phi),

    // <--- Bulirsch's Integrals --->

    /// Compute Bulirsch's complete elliptic integral.
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert. Note that the original
    /// literature by Bulirsch used the complementary modulus kc where kc = √(1 - m).
    ///
    /// ## Parameters
    ///
    /// - `self`: elliptic parameter (m)
    /// - `p`: a parameter
    /// - `a`, `b`: coefficients
    ///
    /// ## Example
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let m = 0.5;
    /// let p = 1.0;
    /// let a = 1.0;
    /// let b = 1.0;
    ///
    /// assert::close(m.cel(p, a, b), 1.8540746773013717, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if m = 1, p = 0, or more than one arguments are infinite.
    ///
    /// [1]: https://crates.io/crates/ellip
    @kc cel -> cel(p, a, b),

    /// Compute Bulirsch's complete elliptic integral.
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert. Note that the original
    /// literature by Bulirsch used the complementary modulus kc where kc = √(1 - m).
    ///
    /// ## Parameters
    ///
    /// - `self`: elliptic parameter (m)
    ///
    /// ## Example
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let m = 0.5;
    ///
    /// assert::close(m.cel1(), 1.8540746773013717, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if m = 1.
    ///
    /// [1]: https://crates.io/crates/ellip
    @kc cel1 -> cel1(),

    /// Compute Bulirsch's complete elliptic integral.
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert. Note that the original
    /// literature by Bulirsch used the complementary modulus kc where kc = √(1 - m).
    ///
    /// ## Parameters
    ///
    /// - `self`: elliptic parameter (m)
    /// - `a`, `b`: coefficients
    ///
    /// ## Example
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let m = 0.5;
    /// let a = 1.0;
    /// let b = 1.0;
    ///
    /// assert::close(m.cel2(a, b), 1.8540746773013717, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if m = 1 or more than one arguments are infinite.
    ///
    /// [1]: https://crates.io/crates/ellip
    @kc cel2 -> cel2(a, b),

    /// Compute Bulirsch's incomplete elliptic integral.
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert. Note that the original
    /// literature by Bulirsch used the complementary modulus kc where kc = √(1 - m).
    ///
    /// ## Parameters
    ///
    /// - `self`: elliptic parameter (m)
    /// - `x`: tangent of amplitude angle
    ///
    /// ## Example
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let m = 0.5;
    /// let x = std::f64::consts::FRAC_PI_4.tan();
    ///
    /// assert::close(m.el1(x), 0.826017876249245, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if m = 1.
    ///
    /// [1]: https://crates.io/crates/ellip
    @kc_x el1 -> el1(x),

    /// Compute Bulirsch's incomplete elliptic integral.
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert. Note that the original
    /// literature by Bulirsch used the complementary modulus kc where kc = √(1 - m).
    ///
    /// ## Parameters
    ///
    /// - `self`: elliptic parameter (m)
    /// - `x`: tangent of amplitude angle
    /// - `a`, `b`: coefficients
    ///
    /// ## Example
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let m = 0.5;
    /// let x = std::f64::consts::FRAC_PI_4.tan();
    /// let a = 1.0;
    /// let b = 1.0;
    ///
    /// assert::close(m.el2(x, a, b), 0.826017876249245, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if m = 1.
    ///
    /// [1]: https://crates.io/crates/ellip
    @kc_x el2 -> el2(x, a, b),

    /// Compute Bulirsch's incomplete elliptic integral.
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert. Note that the original
    /// literature by Bulirsch used the complementary modulus kc where kc = √(1 - m).
    ///
    /// ## Parameters
    ///
    /// - `self`: elliptic parameter (m)
    /// - `x`: tangent of amplitude angle
    /// - `p`: a parameter
    ///
    /// ## Example
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let m = 0.5;
    /// let x = std::f64::consts::FRAC_PI_4.tan();
    /// let p = 1.0;
    ///
    /// assert::close(m.el3(x, p), 0.826017876249245, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if m = 1, 1 + px² = 0, or m < 0 for p < 0.
    ///
    /// [1]: https://crates.io/crates/ellip
    @kc_x el3 -> el3(x, p),

    // <--- Carlson's Symmetric Integrals --->

    /// Compute Carlson's symmetric elliptic integral of the first kind (RF).
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert.
    ///
    /// ## Parameters
    ///
    /// - `self`, `y`, `z`: symmetric arguments
    ///
    /// ## Example
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let x = 1.0;
    /// let y = 0.5;
    /// let z = 0.25;
    ///
    /// assert::close(x.elliprf(y, z), 1.370171633266872, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if any of x, y, or z is negative, or more than one of them are zero.
    ///
    /// [1]: https://crates.io/crates/ellip
    elliprf -> elliprf(y, z),

    /// Compute Carlson's symmetric elliptic integral of the second kind (RG).
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert.
    ///
    /// ## Parameters
    ///
    /// - `self`, `y`, `z`: symmetric arguments
    ///
    /// ## Example
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let x = 1.0;
    /// let y = 0.5;
    /// let z = 0.25;
    ///
    /// assert::close(x.elliprg(y, z), 0.7526721491833781, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if any of x, y, or z is negative or infinite.
    ///
    /// [1]: https://crates.io/crates/ellip
    elliprg -> elliprg(y, z),

    /// Compute Carlson's symmetric elliptic integral of the third kind (RJ).
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert.
    ///
    /// ## Parameters
    ///
    /// - `self`, `y`, `z`: symmetric arguments
    /// - `p`: a parameter
    ///
    /// ## Example
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let x = 1.0;
    /// let y = 0.5;
    /// let z = 0.25;
    /// let p = 0.125;
    ///
    /// assert::close(x.elliprj(y, z, p), 5.680557292035963, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if p = 0, any of x, y, or z is negative, or more than one of them are zero.
    ///
    /// [1]: https://crates.io/crates/ellip
    @first elliprj -> elliprj(y, z, p),

    /// Compute Carlson's degenerate symmetric elliptic integral (RC).
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert.
    ///
    /// ## Parameters
    ///
    /// - `self`, `y`: symmetric arguments
    ///
    /// ## Example
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let x = 1.0;
    /// let y = 0.5;
    ///
    /// assert::close(x.elliprc(y), 1.2464504802804608, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if x < 0, y = 0, or y < 0.
    ///
    /// [1]: https://crates.io/crates/ellip
    @first elliprc -> elliprc(y),

    /// Compute Carlson's symmetric elliptic integral of the second kind (RD).
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert.
    ///
    /// ## Parameters
    ///
    /// - `self`, `y`, `z`: symmetric arguments
    ///
    /// ## Example
    ///
    /// ```
    /// use special::Elliptic;
    ///
    /// let x = 1.0;
    /// let y = 0.5;
    /// let z = 0.25;
    ///
    /// assert::close(x.elliprd(y, z), 4.022594757168912, 1e-15)
    /// ```
    ///
    /// ## Panics
    ///
    /// The function panics if x < 0, y < 0, z ≤ 0 or when both x and y are zero.
    ///
    /// [1]: https://crates.io/crates/ellip
    @first elliprd -> elliprd(y, z),

    // <--- Miscellaneous Functions --->

    /// Compute Jacobi's Zeta function.
    ///
    /// The implementation is based on [ellip][1] by Sira Pornsiriprasert.
    ///
    /// ## Parameters
    ///
    /// - `self`: elliptic parameter (m)
    /// - `phi`: amplitude angle (φ)
    ///
    /// ## Example
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
    /// The function panics if m > 1, phi is infinite, or m is infinite.
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
    /// ## Example
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
