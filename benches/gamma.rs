use random::{self, Source};
use special;
use test::{Bencher, black_box};

#[bench]
fn digamma(b: &mut Bencher) {
    let x = random::default().iter::<f64>().take(1000).map(|x| 20.0 * x - 10.0)
                                                      .collect::<Vec<_>>();

    b.iter(|| {
        for &x in x.iter() {
            black_box(special::digamma(x));
        }
    });
}
