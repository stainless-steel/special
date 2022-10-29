#![feature(test)]

extern crate random;
extern crate special;
extern crate test;

use random::Source;
use special::Gamma;
use test::{black_box, Bencher};

#[bench]
fn digamma(bencher: &mut Bencher) {
    let x = random::default()
        .iter::<f64>()
        .take(1000)
        .map(|x| 20.0 * x - 10.0)
        .collect::<Vec<_>>();
    bencher.iter(|| {
        for &x in &x {
            black_box(x.digamma());
        }
    });
}

#[bench]
fn inc_gamma(bencher: &mut Bencher) {
    let (mut x, mut p) = (random::default(), random::default());
    let xp = x
        .iter::<f64>()
        .zip(p.iter::<f64>())
        .take(1000)
        .map(|(x, p)| (100.0 * x, 100.0 * p))
        .collect::<Vec<_>>();
    bencher.iter(|| {
        for &(x, p) in &xp {
            black_box(x.inc_gamma(p));
        }
    });
}

#[bench]
fn trigamma(bencher: &mut Bencher) {
    let x = random::default().iter().take(1000).collect::<Vec<f64>>();

    bencher.iter(|| {
        for &x in &x {
            black_box(x.trigamma());
        }
    });
}

#[bench]
fn gamma(bencher: &mut Bencher) {
    let rng = random::Default::new();
    let x = rng
        .seed([202208261, 2118])
        .iter::<f64>()
        .take(1000)
        .map(|x| 20.0 * x - 10.0)
        .collect::<Vec<_>>();

    bencher.iter(|| {
        let mut checksum = 0; // See comment in `error.rs`.
        for &x in &x {
            let value = x.gamma();
            checksum ^= unsafe { core::mem::transmute::<f64, u64>(value) };
        }
        black_box(checksum);
    });
}

#[bench]
fn ln_gamma(bencher: &mut Bencher) {
    let rng = random::Default::new();
    let x = rng
        .seed([202208262, 2119])
        .iter::<f64>()
        .take(1000)
        .map(|x| 20.0 * x - 10.0)
        .collect::<Vec<_>>();

    bencher.iter(|| {
        let mut checksum = 0; // See comment in `error.rs`.
        for &x in &x {
            let (value, _sign) = x.ln_gamma();
            checksum ^= unsafe { core::mem::transmute::<f64, u64>(value) };
        }
        black_box(checksum);
    });
}
