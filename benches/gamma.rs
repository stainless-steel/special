#![feature(test)]

extern crate random;
extern crate special;
extern crate test;

use random::Source;
use special::Gamma;
use test::{black_box, Bencher};

#[bench]
fn digamma(bencher: &mut Bencher) {
    let x = random::default(42)
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
    let (mut x, mut p) = (random::default(42), random::default(24));
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
    let x = random::default(42).iter().take(1000).collect::<Vec<f64>>();

    bencher.iter(|| {
        for &x in &x {
            black_box(x.trigamma());
        }
    });
}
