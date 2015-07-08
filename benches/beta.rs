use random::{self, Source};
use special;
use test::{Bencher, black_box};

#[bench]
fn inc_beta(bencher: &mut Bencher) {
    let (p, q) = (0.5, 1.5);
    let ln_beta = special::ln_beta(p, q);
    let x = random::default().iter().take(1000).collect::<Vec<_>>();

    bencher.iter(|| {
        for &x in &x {
            black_box(special::inc_beta(x, p, q, ln_beta));
        }
    });
}

#[bench]
fn inv_inc_beta(bencher: &mut Bencher) {
    let (p, q) = (0.5, 1.5);
    let ln_beta = special::ln_beta(p, q);
    let x = random::default().iter().take(1000).collect::<Vec<_>>();

    bencher.iter(|| {
        for &x in &x {
            black_box(special::inv_inc_beta(x, p, q, ln_beta));
        }
    });
}
