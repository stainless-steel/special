/// Re-export float-point constants defined in the standard library as associated constants.
pub(crate) trait FpConsts {
    const PI: Self;
}

impl FpConsts for f32 {
    const PI: Self = std::f32::consts::PI;
}

impl FpConsts for f64 {
    const PI: Self = std::f64::consts::PI;
}
