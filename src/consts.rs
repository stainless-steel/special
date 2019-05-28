pub(crate) trait Float {
    const PI: Self;
}

impl Float for f32 {
    const PI: Self = std::f32::consts::PI;
}

impl Float for f64 {
    const PI: Self = std::f64::consts::PI;
}
