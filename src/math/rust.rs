#[allow(non_camel_case_types)]
type c_double = f64;

#[allow(non_camel_case_types)]
type c_int = i32;

#[inline(always)]
pub unsafe fn erf(x: c_double) -> c_double {
    libm::erf(x)
}

#[inline(always)]
pub unsafe fn erfc(x: c_double) -> c_double {
    libm::erfc(x)
}

#[inline(always)]
pub unsafe fn tgamma(x: c_double) -> c_double {
    libm::tgamma(x)
}

#[inline(always)]
pub unsafe fn lgamma(x: c_double, sign: &mut c_int) -> c_double {
    let (result, sign_out) = libm::lgamma_r(x);
    *sign = sign_out;
    result
}
