use pairing::bn256::Fq as Fp;
use pairing::group::ff::Field;

// isogeny structure
#[derive(Clone, Debug)]
pub(crate) struct Isogeny {
    a: Fp,
    b: Fp,
}
impl Isogeny {
    pub(crate) fn new(depth: usize) -> Isogeny {
        match depth {
            0 => Isogeny {
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
            1 => Isogeny {
                a: Fp::from_raw([
                    0xff34c82326b2aed2,
                    0x8ff011803978d379,
                    0x961e8c0d4c7d2e60,
                    0x3061dbb2a0530bbc,
                ]),
                b: Fp::from_raw([
                    0x4ba6ec7c7c884f9e,
                    0xdc06b468c744c7a2,
                    0xc0e12fbb87e96dbc,
                    0x218b98fa28f305d7,
                ]),
            },
            2 => Isogeny {
                a: Fp::from_raw([
                    0x79d787e763949cea,
                    0x0f22b2225df1ee26,
                    0x446373526a0853fa,
                    0x4e58081bd28da,
                ]),
                b: Fp::from_raw([
                    0x1717c94a5d56e544,
                    0x59a5b8447109047a,
                    0x95b048d593794895,
                    0x2bfa1dbfff232478,
                ]),
            },
            3 => Isogeny {
                a: Fp::from_raw([
                    0x4269f23016c49098,
                    0x136fac5584322e21,
                    0x93b872de5cefef45,
                    0x15eec1ab01cc5739,
                ]),
                b: Fp::from_raw([
                    0x11c8bfd2ef45ca6f,
                    0x2d2634e2be50ae90,
                    0xd8af4283874212df,
                    0x2d78214c09f651bb,
                ]),
            },
            4 => Isogeny {
                a: Fp::from_raw([
                    0x27c44a827d551f25,
                    0x37d0dd8558faf0c4,
                    0x1142f3515342e9fa,
                    0x24c97155663f5abe,
                ]),
                b: Fp::from_raw([
                    0x8d4c26e8dc988069,
                    0x425e2f3ac3a7b295,
                    0x91a671f76214bad5,
                    0x10cfd9bf7cf2a732,
                ]),
            },
            5 => Isogeny {
                a: Fp::from_raw([
                    0x8521c70ea1542081,
                    0x2c31a2bad4469515,
                    0xe746cb181c953606,
                    0xafd270865b38f5f,
                ]),
                b: Fp::from_raw([
                    0xa1210718a65f2a66,
                    0x8f95ad6dc69793e1,
                    0xd987e468fdc617e2,
                    0x2da42b9297c73294,
                ]),
            },
            6 => Isogeny {
                a: Fp::from_raw([
                    0xcffcdc3142bb4432,
                    0x2fe9c0dfa825dcc3,
                    0xe829d90d0e49eca7,
                    0x1b3aa11377439580,
                ]),
                b: Fp::from_raw([
                    0x6581a03b3c7d5bb0,
                    0x9e6dd048a6e74600,
                    0x394c3a48ba8281d9,
                    0x289ef4b2b9d1a3cd,
                ]),
            },
            7 => Isogeny {
                a: Fp::from_raw([
                    0xad1a5ed6e5b28a8a,
                    0x8a901590dbb2a28a,
                    0xb094b892b928837c,
                    0x15531d00ea03d849,
                ]),
                b: Fp::from_raw([
                    0xcf8a39f73ec5018a,
                    0xcac03fafe3730d68,
                    0x8e2b0180e6af08ba,
                    0x1b411ddaf2fea5a7,
                ]),
            },
            8 => Isogeny {
                a: Fp::from_raw([
                    0x08ad0e3aff97ba9b,
                    0xf73d3975ac88cc05,
                    0xbff7c0642e5f573f,
                    0x219304605e89f38e,
                ]),
                b: Fp::from_raw([
                    0xef3fff59344445a2,
                    0x8684c6b5dc31b97f,
                    0x97ab0877c200a627,
                    0x2138d0df886dc840,
                ]),
            },
            9 => Isogeny {
                a: Fp::from_raw([
                    0xb5e2c75f6ce2d930,
                    0xafe1399906152298,
                    0x70e75fcfbd5a1b92,
                    0x3fb95fe2f9dccd7,
                ]),
                b: Fp::from_raw([
                    0x4b48eede3c8f2b5a,
                    0x96ee646f699f7d9d,
                    0xfe538433e387a6b2,
                    0x2a3ca87641d9e5b2,
                ]),
            },
            10 => Isogeny {
                a: Fp::from_raw([
                    0x484137768937949b,
                    0x5e2c9ce6a7cff768,
                    0xcd61bd13d773e28b,
                    0x1dfb93ae3ff8048,
                ]),
                b: Fp::from_raw([
                    0xac3ca96707aad4ea,
                    0x82c499b8547437da,
                    0xa9f51753e4f8afee,
                    0x13f681e31e6cfd15,
                ]),
            },
            11 => Isogeny {
                a: Fp::from_raw([
                    0x438ba164c5dfc246,
                    0xd801a0436928dce4,
                    0xae3818b317ed6d81,
                    0xf47a855cb7e2c1,
                ]),
                b: Fp::from_raw([
                    0x3159473582891908,
                    0x252d7e8b9f413cb5,
                    0x7d3e6ab3ca2e61b0,
                    0x226282adc850179b,
                ]),
            },
            12 => Isogeny {
                a: Fp::from_raw([
                    0xa2b052f952b6a765,
                    0xf79fb9d58df8ca3b,
                    0xfc9d439158d4b4f2,
                    0x1410f01c307d3654,
                ]),
                b: Fp::from_raw([
                    0xe91de3ee372d4ae3,
                    0x296ba83eec84f8a5,
                    0x4356b81b6e192a08,
                    0x1665b1b81ae65291,
                ]),
            },
            13 => Isogeny {
                a: Fp::from_raw([
                    0x1bc20bf2f984b890,
                    0xaa21391d073eec4d,
                    0xc0f76b4d112c04df,
                    0xbe54c10dbee3311,
                ]),
                b: Fp::from_raw([
                    0x93623bd80a78695a,
                    0x4fa46fa95dbec430,
                    0xf06c6d97728fe37a,
                    0x2a0aa2438baedb87,
                ]),
            },
            _ => {
                unimplemented!()
            }
        }
    }

    pub(crate) fn evaluate(&self, x: Fp) -> Fp {
        let Isogeny { a, b } = self;
        let numerator = *b + *a * x + x.square();
        let denominator = *a + x;
        numerator * denominator.invert().unwrap()
    }
}
