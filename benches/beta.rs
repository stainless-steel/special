use random::{self, Source};
use special::Beta;
use test::{Bencher, black_box};

#[bench]
fn inc_beta(bencher: &mut Bencher) {
    let (p, q) = (0.5, 1.5);
    let ln_beta = p.ln_beta(q);
    let x = random::default().iter().take(1000).collect::<Vec<f64>>();

    bencher.iter(|| {
        for &x in &x {
            black_box(x.inc_beta(p, q, ln_beta));
        }
    });
}

#[bench]
fn inv_inc_beta(bencher: &mut Bencher) {
    let (p, q) = (0.5, 1.5);
    let ln_beta = p.ln_beta(q);
    let alpha = random::default().iter().take(1000).collect::<Vec<f64>>();

    bencher.iter(|| {
        for &alpha in &alpha {
            black_box(alpha.inv_inc_beta(p, q, ln_beta));
        }
    });
}
