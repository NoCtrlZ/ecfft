use pairing::arithmetic::BaseExt;
use pairing::bn256::Fq as Fp;
use pairing::group::ff::Field;
use rayon::prelude::*;

// isogeny structure
#[derive(Clone, Debug)]
pub(crate) struct Isogeny {
    a: Fp,
    b: Fp,
}
impl Isogeny {
    pub(crate) fn new(depth: usize) -> Isogeny {
        match depth {
            1 => Isogeny {
                a: Fp::from_raw([
                    0x3c8628054523a5de,
                    0x4f8961d14bb560d0,
                    0x6d40ffafdb42c12d,
                    0x1833609991081a4b,
                ]),
                b: Fp::from_raw([
                    0xf0799f9a5bf4ada6,
                    0xbb7ab628a12d02ea,
                    0xf76f15faf997eaa0,
                    0xed8b578b83e9a51,
                ]),
            },
            2 => Isogeny {
                a: Fp::from_raw([
                    0x509a7c8c05b12426,
                    0x44dbeb15610c8b88,
                    0x64ee1cb7973bfbd1,
                    0x57BB06ac07315ce,
                ]),
                b: Fp::from_raw([
                    0xcb6ac9472dab0b76,
                    0xe51b01edc996d386,
                    0xfe2e12a8111cb7d6,
                    0x1803646703251b2d,
                ]),
            },
            3 => Isogeny {
                a: Fp::from_raw([
                    0x58f935a65574871b,
                    0xf3d49205b05b2ed4,
                    0xf264ce41f5311095,
                    0x154b6ff211dc3eb9,
                ]),
                b: Fp::from_raw([
                    0x2c96cb2ffb5157db,
                    0x2d9df99cc2c197d2,
                    0xc49f6b7ade396133,
                    0x413428325e24475,
                ]),
            },
            4 => Isogeny {
                a: Fp::from_raw([
                    0xb05094c95e744772,
                    0xf0ecc3530f2e17e8,
                    0x67e5c433a785a398,
                    0xed85d5ed1b94be2,
                ]),
                b: Fp::from_raw([
                    0x109e44fa1b94d1a1,
                    0x41c9e2cd63942533,
                    0xf2b6986b006ce2a1,
                    0x14ffe0243def0f38,
                ]),
            },
            5 => Isogeny {
                a: Fp::from_raw([
                    0x920f7d17bced4fb0,
                    0x97bb25809e425c77,
                    0x1632991e84532758,
                    0x1f00cf7e4e69b575,
                ]),
                b: Fp::from_raw([
                    0x06581a03b3c7d5bb,
                    0x99e6dd048a6e7460,
                    0xd394c3a48ba8281d,
                    0x289ef4b2b9d1a3c,
                ]),
            },
            6 => Isogeny {
                a: Fp::from_raw([
                    0x4956ddc125ab2146,
                    0xee64baaceb258de9,
                    0x484d50ffef0acd0d,
                    0x1d86ee79ab19c627,
                ]),
                b: Fp::from_raw([
                    0xd28cfb2dbb3a6e65,
                    0x2b5ce695df7e4f6e,
                    0x8c14dbaa1f5bc7c6,
                    0x1ff2c2e57beeee74,
                ]),
            },
            7 => Isogeny {
                a: Fp::from_raw([
                    0x2f43ac9fe243ac9c,
                    0x2f705e4a79778aeb,
                    0xfa3a2461ecb8d816,
                    0x2caffbee4087b502,
                ]),
                b: Fp::from_raw([
                    0xa67811786e53e403,
                    0x2b5879bd8ad154e9,
                    0x4084b93e4c50356e,
                    0x82016dc54ad1089,
                ]),
            },
            8 => Isogeny {
                a: Fp::from_raw([
                    0x2d78b1d7db38b64c,
                    0xabf84e66418548a6,
                    0xdc39d7f3ef5686e4,
                    0xfee57f8be77335,
                ]),
                b: Fp::from_raw([
                    0x1a48e67c2b171102,
                    0xa81fc8e1d7e11672,
                    0x431763d54f2951a5,
                    0x20e27b8f30dca275,
                ]),
            },
            9 => Isogeny {
                a: Fp::from_raw([
                    0xff28b6eec4aba31c,
                    0x092c3726b84955c3,
                    0x7d94a38dd6fdfae9,
                    0x24c32924e1e51831,
                ]),
                b: Fp::from_raw([
                    0xf0582224b7c8cb9b,
                    0x66dd2c36668e6215,
                    0x6dd17d074f406239,
                    0x1f7e1925fea5d3eb,
                ]),
            },
            10 => Isogeny {
                a: Fp::from_raw([
                    0xaef32e649db66f35,
                    0x41c11d598e831c7f,
                    0x47b6290806bc078f,
                    0x186f45dac7c6c8c5,
                ]),
                b: Fp::from_raw([
                    0xf125da7ec4671034,
                    0xce138d316e2cf911,
                    0x93fc09867d639249,
                    0x1a584f644d1dd18e,
                ]),
            },
            11 => Isogeny {
                a: Fp::from_raw([
                    0x37b437c40acce92b,
                    0x23c84919bd9aa532,
                    0xad3b6251f6958354,
                    0x111d4fa3c46bb59f,
                ]),
                b: Fp::from_raw([
                    0x77e83e8e9848c2cf,
                    0xdabfb3c7e6968acb,
                    0xc2ec9b6f2fea7f60,
                    0x22ab510a7c808345,
                ]),
            },
            12 => Isogeny {
                a: Fp::from_raw([
                    0x46f082fcbe612e24,
                    0xea884e4741cfbb13,
                    0x703ddad3444b0137,
                    0x2f9530436fb8cc4,
                ]),
                b: Fp::from_raw([
                    0x4eca7b4bc7f5a4e2,
                    0x23ab299577230adb,
                    0x9238f26b8819d572,
                    0x20df5b2c0579f1d2,
                ]),
            },
            13 => Isogeny {
                a: Fp::from_raw([
                    0xbc0bcabefcf1af9c,
                    0x7ba5bdbe4ee0676b,
                    0x8011a114208a4c36,
                    0x2ab8998f464eee19,
                ]),
                b: Fp::from_raw([
                    0x823df7ecfb045196,
                    0x45be4532c6d517ae,
                    0xc9dcdde38996b870,
                    0x9813e37112d1cf3,
                ]),
            },
            _ => Isogeny {
                a: Fp::one(),
                b: Fp::one(),
            },
        }
    }

    pub(crate) fn evaluate(&self, x: Fp) -> Fp {
        let Isogeny { a, b } = self;
        let numerator = *b + *a * x + x.square();
        let denominator = *a + x;
        numerator * denominator.invert().unwrap()
    }

    pub(crate) fn evaluate_with_denominator(&self, x: Fp) -> Fp {
        let Isogeny { a, b: _ } = self;
        *a + x
    }

    pub(crate) fn domain_half_sizing(&self, domain: Vec<Fp>, size: usize) -> Vec<Fp> {
        domain[..size]
            .par_iter()
            .map(|coeff| self.evaluate(*coeff))
            .collect()
    }

    pub(crate) fn get_factor(
        &self,
        domain: &Vec<Fp>,
        size: usize,
        exp: &[u64; 4],
    ) -> Vec<((Fp, Fp), (Fp, Fp))> {
        domain[..size]
            .par_iter()
            .zip(&domain[size..])
            .map(|(a, b)| {
                let f1 = self.evaluate_with_denominator(*a).pow(exp);
                let f2 = a * f1;
                let f3 = self.evaluate_with_denominator(*b).pow(exp);
                let f4 = b * f3;
                ((f1, f2), (f3, f4))
            })
            .collect()
    }
}
