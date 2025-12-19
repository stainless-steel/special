/// Elliptic integral functions
macro_rules! interface_method {
    ($(#[$attr:meta])* $fn_name:ident $($params:ident)*) => {
        $(#[$attr])*
        fn $fn_name(self, $($params: Self,)*) -> Self;
    };
}

macro_rules! m_to_kc {
    ($fn_name:ident, $m:expr) => {{
        if $m < 0.0 {
            panic!(concat!(
                stringify!($fn_name),
                ": m (self) cannot be less than 1."
            ))
        }
        (1.0 - $m).sqrt()
    }};
}

// By default, impl_method expect self as the last argument.
// With @kc, impl_method expect self as **m** as the first argument.
// With @kc_x, impl_method expect x as the first and self as the second argument.
// With @first, impl_method expect self as the first argument.
macro_rules! impl_method {
    (@kc $fn_name:ident $($params:ident)*) => {
        fn $fn_name(self, $($params: Self,)*) -> Self {
            let kc = m_to_kc!($fn_name, self);
            ellip::$fn_name(kc, $($params,)*).unwrap()
        }
    };
    (@kc_x $fn_name:ident $x:ident) => {
        fn $fn_name(self, $x: Self) -> Self {
            let kc = m_to_kc!($fn_name, self);
            ellip::$fn_name($x, kc).unwrap()
        }
    };
    (@kc_x $fn_name:ident $x:ident $($params:ident)+) => {
        fn $fn_name(self, $x: Self, $($params: Self,)*) -> Self {
            let kc = m_to_kc!($fn_name, self);
            ellip::$fn_name($x, kc, $($params,)*).unwrap()
        }
    };
    (@first $fn_name:ident $($params:ident)*) => {
        fn $fn_name(self, $($params: Self,)*) -> Self {
            ellip::$fn_name(self, $($params,)*).unwrap()
        }
    };
    ($fn_name:ident $($params:ident)*) => {
        fn $fn_name(self, $($params: Self,)*) -> Self {
            ellip::$fn_name($($params,)* self).unwrap()
        }
    };
}

macro_rules! impl_elliptic {
    ($(
        $(comment:meta)*
        $(#[$attr:meta])*
        $(@$flag:ident)? $fn_name:ident $([$($order:tt)*])? $($params:ident)*,
    )*) => {
        /// Elliptic integrals.
        pub trait Elliptic: Sized {
            $(
                interface_method!($(#[$attr])* $fn_name $($params)*);
            )*
        }

        impl Elliptic for f32 {
            $(
                impl_method!($(@$flag)? $fn_name $($params)* );
            )*
        }

        impl Elliptic for f64 {
            $(
                impl_method!($(@$flag)? $fn_name $($params)* );
            )*
        }
    };
}

impl_elliptic!(
    // <--- Complete Legendre's Integrals --->

    /// Computes complete elliptic integral of the first kind (K).
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
    /// ## Domain
    ///
    /// - Returns error if m > 1.
    ///
    /// Implementation based on [Ellip](https://crates.io/crates/ellip) by Sira Pornsiriprasert.
    ellipk,

    /// Computes complete elliptic integral of the second kind (E).
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
    /// ## Domain
    ///
    /// - Returns error if m > 1.
    ///
    /// Implementation based on [Ellip](https://crates.io/crates/ellip) by Sira Pornsiriprasert.
    ellipe,

    /// Computes complete elliptic integral of the third kind (Π).
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
    /// ## Domain
    ///
    /// - Returns error if m > 1.
    /// - Returns the Cauchy principal value if n > 1.
    ///
    /// Implementation based on [Ellip](https://crates.io/crates/ellip) by Sira Pornsiriprasert.
    ellippi n,

    /// Computes complete elliptic integral of Legendre's type (D).
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
    /// ## Domain
    ///
    /// - Returns error if m > 1.
    ///
    /// Implementation based on [Ellip](https://crates.io/crates/ellip) by Sira Pornsiriprasert.
    ellipd,

    // <--- Incomplete Legendre's Integrals --->

    /// Computes incomplete elliptic integral of the first kind (F).
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
    /// ## Domain
    ///
    /// - Returns error if m sin²φ > 1.
    ///
    /// Implementation based on [Ellip](https://crates.io/crates/ellip) by Sira Pornsiriprasert.
    ellipf phi,

    /// Computes incomplete elliptic integral of the second kind (E).
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
    /// ## Domain
    ///
    /// - Returns error if m sin²φ > 1.
    ///
    /// Implementation based on [Ellip](https://crates.io/crates/ellip) by Sira Pornsiriprasert.
    ellipeinc phi,

    /// Computes incomplete elliptic integral of the third kind (Π).
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
    /// ## Domain
    ///
    /// - Returns error when:
    ///   - m sin²φ > 1,
    ///   - n sin²φ = 1,
    ///   - or m ≥ 1 and φ is not a multiple of π/2.
    /// - Returns the Cauchy principal value if n sin²φ > 1.
    ///
    /// Implementation based on [Ellip](https://crates.io/crates/ellip) by Sira Pornsiriprasert.
    ellippiinc phi n,

    /// Computes incomplete elliptic integral of the third kind (Π) using Bulirsch's method.
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
    /// ## Domain
    ///
    /// - Returns error when:
    ///   - m sin²φ > 1,
    ///   - n sin²φ = 1,
    ///   - or m ≥ 1 and φ is not a multiple of π/2.
    /// - Returns the Cauchy principal value if n sin²φ > 1.
    ///
    /// Implementation based on [Ellip](https://crates.io/crates/ellip) by Sira Pornsiriprasert.
    ellippiinc_bulirsch phi n,

    /// Computes incomplete elliptic integral of Legendre's type (D).
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
    /// ## Domain
    ///
    /// - Returns error if m sin²φ > 1.
    ///
    /// Implementation based on [Ellip](https://crates.io/crates/ellip) by Sira Pornsiriprasert.
    ellipdinc phi,

    // <--- Bulirsch's Integrals --->

    /// Computes Bulirsch's complete elliptic integral.
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
    /// ## Domain
    ///
    /// - Returns error if m = 1 or p = 0.
    /// - Returns the Cauchy principal value for p < 0.
    /// - Returns error if more than one arguments are infinite.
    ///
    /// ## Notes
    ///
    /// The original literature by Bulirsch used the complementary modulus `kc`, where kc = √(1-m).
    ///
    /// Implementation based on [Ellip](https://crates.io/crates/ellip) by Sira Pornsiriprasert.
    @kc cel p a b,

    /// Computes Bulirsch's complete elliptic integral.
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
    /// ## Domain
    ///
    /// - Returns error if m = 1.
    ///
    /// ## Notes
    ///
    /// The original literature by Bulirsch used the complementary modulus `kc`, where kc = √(1-m).
    ///
    /// Implementation based on [Ellip](https://crates.io/crates/ellip) by Sira Pornsiriprasert.
    @kc cel1,

    /// Computes Bulirsch's complete elliptic integral.
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
    /// ## Domain
    ///
    /// - Returns error if m = 1.
    /// - Returns error if more than one arguments are infinite.
    ///
    /// ## Notes
    ///
    /// The original literature by Bulirsch used the complementary modulus `kc`, where kc = √(1-m).
    ///
    /// Implementation based on [Ellip](https://crates.io/crates/ellip) by Sira Pornsiriprasert.
    @kc cel2 a b,

    /// Computes Bulirsch's incomplete elliptic integral.
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
    /// ## Domain
    ///
    /// - Returns error if m = 1.
    ///
    /// ## Notes
    ///
    /// The original literature by Bulirsch used the complementary modulus `kc`, where kc = √(1-m).
    ///
    /// Implementation based on [Ellip](https://crates.io/crates/ellip) by Sira Pornsiriprasert.
    @kc_x el1 x,

    /// Computes Bulirsch's incomplete elliptic integral.
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
    /// ## Domain
    ///
    /// - Returns error if m = 1.
    ///
    /// ## Notes
    ///
    /// The original literature by Bulirsch used the complementary modulus `kc`, where kc = √(1-m).
    ///
    /// Implementation based on [Ellip](https://crates.io/crates/ellip) by Sira Pornsiriprasert.
    @kc_x el2 x a b,

    /// Computes Bulirsch's incomplete elliptic integral.
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
    /// ## Domain
    ///
    /// - Returns error if:
    ///   - m = 1,
    ///   - 1 + px² = 0,
    ///   - or m < 0 for p < 0.
    /// - Returns the Cauchy principal value when 1 + px² < 0
    ///
    /// ## Notes
    ///
    /// The original literature by Bulirsch used the complementary modulus `kc`, where kc = √(1-m).
    ///
    /// Implementation based on [Ellip](https://crates.io/crates/ellip) by Sira Pornsiriprasert.
    @kc_x el3 x p,

    // <--- Carlson's Symmetric Integrals --->

    /// Computes Carlson's symmetric elliptic integral of the first kind (RF).
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
    /// ## Domain
    ///
    /// - Returns error if any of x, y, or z is negative, or more than one of them are zero.
    ///
    /// Implementation based on [Ellip](https://crates.io/crates/ellip) by Sira Pornsiriprasert.
    elliprf y z,

    /// Computes Carlson's symmetric elliptic integral of the second kind (RG).
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
    /// ## Domain
    ///
    /// - Returns error if any of x, y, or z is negative or infinite.
    ///
    /// Implementation based on [Ellip](https://crates.io/crates/ellip) by Sira Pornsiriprasert.
    elliprg y z,

    /// Computes Carlson's symmetric elliptic integral of the third kind (RJ).
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
    /// ## Domain
    ///
    /// - Returns error if:
    ///   - any of x, y, or z is negative, or more than one of them are zero,
    ///   - or p = 0.
    /// - Returns the Cauchy principal value if p < 0.
    ///
    /// Implementation based on [Ellip](https://crates.io/crates/ellip) by Sira Pornsiriprasert.
    @first elliprj y z p,

    /// Computes Carlson's degenerate symmetric elliptic integral (RC).
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
    /// ## Domain
    ///
    /// - Returns error if x < 0 or y = 0.
    /// - Returns the Cauchy principal value if y < 0.
    ///
    /// Implementation based on [Ellip](https://crates.io/crates/ellip) by Sira Pornsiriprasert.
    @first elliprc y,

    /// Computes Carlson's symmetric elliptic integral of the second kind (RD).
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
    /// ## Domain
    ///
    /// - Returns error if x < 0, y < 0, z ≤ 0 or when both x and y are zero.
    ///
    /// Implementation based on [Ellip](https://crates.io/crates/ellip) by Sira Pornsiriprasert.
    @first elliprd y z,

    // <--- Miscellaneous Functions --->

    /// Computes Jacobi's Zeta function.
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
    /// ## Domain
    ///
    /// - Returns error if m > 1.
    /// - Returns error if phi or m is infinite.
    ///
    /// Implementation based on [Ellip](https://crates.io/crates/ellip) by Sira Pornsiriprasert.
    jacobi_zeta phi,

    /// Computes Heuman's Lambda function.
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
    /// ## Domain
    ///
    /// - Returns error if m < 0 or m ≥ 1.
    /// - Returns error if phi is infinite.
    ///
    /// Implementation based on [Ellip](https://crates.io/crates/ellip) by Sira Pornsiriprasert.
    heuman_lambda phi,
);
