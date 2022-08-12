use pairing::bn256::Fq as Fp;
use pairing::group::ff::Field;

// isogeny structure
#[derive(Clone, Debug)]
pub struct Isogeny {
    pub numerator: [Fp; 2],
    pub denominator: [Fp; 1],
}

impl Isogeny {
    pub(crate) fn evaluate(&self, x: Fp) -> Fp {
        let Isogeny {
            numerator: [a0, a1],
            denominator: [b0],
        } = self;
        let numerator = *a0 + *a1 * x + x * x;
        let denominator = *b0 + x;
        numerator * denominator.invert().unwrap()
    }
}

pub(crate) fn get_isogeny(depth: usize) -> Isogeny {
    match depth {
        1 => Isogeny {
            numerator: [
                Fp::from_raw([
                    0xff34c82326b2aed2,
                    0x8ff011803978d379,
                    0x961e8c0d4c7d2e60,
                    0x3061dbb2a0530bbc,
                ]),
                Fp::from_raw([
                    0x4ba6ec7c7c884f9e,
                    0xdc06b468c744c7a2,
                    0xc0e12fbb87e96dbc,
                    0x218b98fa28f305d7,
                ]),
            ],
            denominator: [Fp::from_raw([
                0xff34c82326b2aed2,
                0x8ff011803978d379,
                0x961e8c0d4c7d2e60,
                0x3061dbb2a0530bbc,
            ])],
        },
        _ => Isogeny {
            numerator: [Fp::zero(), Fp::zero()],
            denominator: [Fp::zero()],
        },
    }
}
