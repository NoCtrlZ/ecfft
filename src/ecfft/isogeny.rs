use pairing::bn256::Fq as Fp;
use pairing::group::ff::Field;

// isogeny structure
#[derive(Clone, Debug)]
pub(crate) struct Isogeny {
    intermidiate: Fp,
    constant: Fp,
}
impl Isogeny {
    pub(crate) fn new(depth: usize) -> Isogeny {
        match depth {
            0 => Isogeny {
                intermidiate: Fp::from_raw([
                    0x3c8628054523a5de,
                    0x4f8961d14bb560d0,
                    0x6d40ffafdb42c12d,
                    0x1833609991081a4b,
                ]),
                constant: Fp::from_raw([
                    0xf0799f9a5bf4ada6,
                    0xbb7ab628a12d02ea,
                    0xf76f15faf997eaa0,
                    0xed8b578b83e9a51,
                ]),
            },
            1 => Isogeny {
                intermidiate: Fp::from_raw([
                    0xff34c82326b2aed2,
                    0x8ff011803978d379,
                    0x961e8c0d4c7d2e60,
                    0x3061dbb2a0530bbc,
                ]),
                constant: Fp::from_raw([
                    0x4ba6ec7c7c884f9e,
                    0xdc06b468c744c7a2,
                    0xc0e12fbb87e96dbc,
                    0x218b98fa28f305d7,
                ]),
            },
            _ => Isogeny {
                intermidiate: Fp::one(),
                constant: Fp::one(),
            },
        }
    }

    pub(crate) fn evaluate(&self, x: Fp) -> Fp {
        let Isogeny {
            intermidiate: c1,
            constant: c2,
        } = self;
        let numerator = *c2 + *c1 * x + x.square();
        let denominator = *c1 + x;
        numerator * denominator.invert().unwrap()
    }
}
