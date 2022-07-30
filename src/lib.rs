mod classic_fft;

pub use classic_fft::ClassicFft;

#[cfg(test)]
mod tests {
    use super::ClassicFft;
    use ff::Field;
    use pasta_curves::Fp;
    use rand_core::OsRng;

    #[test]
    fn fft_test() {
        let k = 12;
        let rng = OsRng;
        let fft = ClassicFft::new(k);
        let mut poly = (0..1 << k).map(|_| Fp::random(rng)).collect::<Vec<_>>();
        fft.fft(&mut poly)
    }
}
