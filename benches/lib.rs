#![feature(test)]

extern crate random;
extern crate special;
extern crate test;

use random::Source;

#[bench]
fn digamma(b: &mut test::Bencher) {
    let x = random::default().iter::<f64>().take(1000).map(|x| 20.0 * x - 10.0)
                                                      .collect::<Vec<_>>();

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
    let x = random::default().iter().take(1000).collect::<Vec<_>>();

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
    let x = random::default().iter().take(1000).collect::<Vec<_>>();

    b.iter(|| {
        for &x in x.iter() {
            test::black_box(special::inv_inc_beta(x, p, q, ln_beta));
        }
    });
}
