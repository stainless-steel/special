#![feature(test)]

extern crate random;
extern crate special;
extern crate test;

use random::Source;
use special::Error;
use test::{black_box, Bencher};

#[bench]
fn error(bencher: &mut Bencher) {
    let rng = random::Default::new();
    let x = rng
        .seed([202208263, 2120])
        .iter::<f64>()
        .take(1000)
        .map(|x| 20.0 * x - 10.0)
        .collect::<Vec<_>>();

    bencher.iter(|| {
        // We avoid calling `black_box` in an inner loop when measuring execution time
        // because it might interact unnaturally with speculative execution. Instead, we
        // calculate an extremely cheap checksum of all calculations and then call
        // `black_box` on the result.
        let mut checksum = 0;
        for &x in &x {
            let value = x.error();
            checksum ^= unsafe { core::mem::transmute::<f64, u64>(value) };
        }
        black_box(checksum);
    });
}

#[bench]
fn compl_error(bencher: &mut Bencher) {
    let rng = random::Default::new();
    let x = rng
        .seed([202208264, 2121])
        .iter::<f64>()
        .take(1000)
        .map(|x| 20.0 * x - 10.0)
        .collect::<Vec<_>>();

    bencher.iter(|| {
        let mut checksum = 0; // See comment in `error` benchmark
        for &x in &x {
            let value = x.compl_error();
            checksum ^= unsafe { core::mem::transmute::<f64, u64>(value) };
        }
        black_box(checksum);
    });
}
