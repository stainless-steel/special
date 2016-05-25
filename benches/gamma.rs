use random::{self, Source};
use special::{self, Gamma};
use test::{Bencher, black_box};

#[bench]
fn digamma(bencher: &mut Bencher) {
    let x = random::default().iter::<f64>().take(1000).map(|x| 20.0 * x - 10.0)
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
    let xp = x.iter::<f64>().zip(p.iter::<f64>()).take(1000).map(|(x, p)| {
        (100.0 * x, 100.0 * p)
    }).collect::<Vec<_>>();

    bencher.iter(|| {
        for &(x, p) in &xp {
            black_box(special::inc_gamma(x, p));
        }
    });
}
