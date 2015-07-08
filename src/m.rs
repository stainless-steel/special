use libc::{c_double, c_int};

extern {
    pub fn erf(x: c_double) -> c_double;
    pub fn erfc(x: c_double) -> c_double;
    pub fn lgamma_r(x: c_double, sign: &mut c_int) -> c_double;
}
