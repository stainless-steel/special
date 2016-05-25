#![link_name = "m"]

use libc::{c_double, c_int};

extern {
    pub fn erf(x: c_double) -> c_double;
    pub fn erfc(x: c_double) -> c_double;
    pub fn tgamma(x: c_double) -> c_double;
}

#[cfg(windows)]
extern {
    pub fn lgamma(x: c_double, sign: &mut c_int) -> c_double;
}

#[cfg(unix)]
extern {
    pub fn lgamma_r(x: c_double, sign: &mut c_int) -> c_double;
}

#[cfg(unix)]
#[inline(always)]
pub unsafe fn lgamma(x: c_double, sign: &mut c_int) -> c_double {
    lgamma_r(x, sign)
}
