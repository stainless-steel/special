#![feature(test)]

extern crate rand;
extern crate test;

extern crate special;

use rand::random;

#[bench]
fn digamma(b: &mut test::Bencher) {
    let x = (0..1000).map(|_| 20.0 * random::<f64>() - 10.0).collect::<Vec<_>>();

    b.iter(|| {
        for &x in x.iter() {
            test::black_box(special::digamma(x));
        }
    });
}

#[bench]
fn inc_beta(b: &mut test::Bencher) {
    let (p, q) = (0.5, 1.5);
    let ln_beta = special::ln_beta(p, q);
    let x = (0..1000).map(|_| random()).collect::<Vec<_>>();

    b.iter(|| {
        for &x in x.iter() {
            test::black_box(special::inc_beta(x, p, q, ln_beta));
        }
    });
}

#[bench]
fn inv_inc_beta(b: &mut test::Bencher) {
    let (p, q) = (0.5, 1.5);
    let ln_beta = special::ln_beta(p, q);
    let x = (0..1000).map(|_| random()).collect::<Vec<_>>();

    b.iter(|| {
        for &x in x.iter() {
            test::black_box(special::inv_inc_beta(x, p, q, ln_beta));
        }
    });
}
