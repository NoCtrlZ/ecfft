use pairing::bn256::Fq as Fp;
use pairing::group::ff::Field;

// isogeny structure
#[derive(Clone, Debug)]
pub(crate) struct Isogeny {
    pub numerator: [Fp; 3],
    pub denominator: [Fp; 2],
}

impl Isogeny {
    pub(crate) fn new(depth: usize) -> Isogeny {
        match depth {
            1 => Isogeny {
                numerator: [
                    Fp::from_raw([
                        0x4ba6ec7c7c884f9e,
                        0xdc06b468c744c7a2,
                        0xc0e12fbb87e96dbc,
                        0x218b98fa28f305d7,
                    ]),
                    Fp::from_raw([
                        0xff34c82326b2aed2,
                        0x8ff011803978d379,
                        0x961e8c0d4c7d2e60,
                        0x3061dbb2a0530bbc,
                    ]),
                    Fp::one(),
                ],
                denominator: [
                    Fp::from_raw([
                        0xff34c82326b2aed2,
                        0x8ff011803978d379,
                        0x961e8c0d4c7d2e60,
                        0x3061dbb2a0530bbc,
                    ]),
                    Fp::one(),
                ],
            },
            _ => Isogeny {
                numerator: [Fp::one(), Fp::one(), Fp::one()],
                denominator: [Fp::one(), Fp::one()],
            },
        }
    }

    pub(crate) fn evaluate(&self, x: Fp) -> Fp {
        let Isogeny {
            numerator: [a0, a1, a2],
            denominator: [b0, b1],
        } = self;
        let numerator = *a0 + *a1 * x + *a2 * x.square();
        let denominator = *b0 + b1 * x;
        numerator * denominator.invert().unwrap()
    }
}
