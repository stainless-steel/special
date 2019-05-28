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

fn trigamma(value: f64) -> f64 {

    const C_LIMIT: f64 = 49.0;
    const S_LIMIT: f64 = 1e-5;

    let mut sum: f64 = 0.0;
    let mut x: f64 = value;

    loop {
        if f64::is_nan(x) || f64::is_infinite(x) {
            return x;
        }

        if x > 0.0 && x <= S_LIMIT {
            return (1.0 / (x * x)) + sum;
        }

        if x >= C_LIMIT {
            let inv = 1.0 / (x * x);
            return (1.0 / x + inv / 2.0 + inv / x * (1.0 / 6.0 - inv * (1.0 / 30.0 + inv / 42.0))) + sum;
        }

        sum += 1.0 / (x * x);
        x += 1.0;
    }
}

#[bench]
fn it_trigamma(bencher: &mut Bencher) {
    let x = random::default().seed([1; 2]).iter().take(1000).collect::<Vec<f64>>();

    bencher.iter(|| {
        for &x in &x {
            black_box(trigamma(x));
        }
    });
}

#[bench]
fn pl_trigamma(bencher: &mut Bencher) {
    let x = random::default().seed([1; 2]).iter().take(1000).collect::<Vec<f64>>();

    bencher.iter(|| {
        for &x in &x {
            black_box(x.trigamma());
        }
    });
}
