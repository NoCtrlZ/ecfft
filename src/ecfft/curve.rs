use core::fmt::{self, Debug};

use ff::Field;
use pasta_curves::Fp;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq, CtOption};

#[derive(Clone)]
pub(crate) struct EpAffine {
    x: Fp,
    y: Fp,
}

impl EpAffine {
    const fn curve_constant_a() -> Fp {
        Fp::from_raw([1, 0, 0, 0])
    }

    const fn curve_constant_b() -> Fp {
        Fp::from_raw([
            0xedc87ab655e55142,
            0xd76d8e4277cb9048,
            0xc6ad51a6a7fe7a43,
            0x34524f71a21a7096,
        ])
    }

    fn generator() -> Self {
        EpAffine {
            x: Fp::from_raw([
                0x4da26202dffa62a8,
                0x8406aac002f8b832,
                0xf60aecbfc30e57f7,
                0x1d62ba1b544c4f84,
            ]),
            y: Fp::from_raw([
                0x52f8ada18c2b96dc,
                0xedd674d51f009506,
                0x17abe167e59849de,
                0x620e16d2e51fdfa,
            ]),
        }
    }

    fn is_on_curve(&self) -> Choice {
        // y^2 - x^3 - ax ?= b
        (self.y.square() - (self.x.square() + Self::curve_constant_a()) * self.x)
            .ct_eq(&Self::curve_constant_b())
            | self.is_identity()
    }

    fn is_identity(&self) -> Choice {
        self.x.is_zero() & self.y.is_zero()
    }
}

impl fmt::Debug for EpAffine {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        if self.is_identity().into() {
            write!(f, "Infinity")
        } else {
            write!(f, "({:?}, {:?})", self.x, self.y)
        }
    }
}

#[derive(Clone, Debug)]
pub(crate) struct Ep {
    x: Fp,
    y: Fp,
    z: Fp,
}

impl Ep {
    const fn curve_constant_a() -> Fp {
        Fp::from_raw([1, 0, 0, 0])
    }

    const fn curve_constant_b() -> Fp {
        Fp::from_raw([
            0xedc87ab655e55142,
            0xd76d8e4277cb9048,
            0xc6ad51a6a7fe7a43,
            0x34524f71a21a7096,
        ])
    }

    fn generator() -> Self {
        Self {
            x: Fp::from_raw([
                0x4da26202dffa62a8,
                0x8406aac002f8b832,
                0xf60aecbfc30e57f7,
                0x1d62ba1b544c4f84,
            ]),
            y: Fp::from_raw([
                0x52f8ada18c2b96dc,
                0xedd674d51f009506,
                0x17abe167e59849de,
                0x620e16d2e51fdfa,
            ]),
            z: Fp::from_raw([1, 0, 0, 0]),
        }
    }

    fn identity() -> Self {
        Self {
            x: Fp::zero(),
            y: Fp::zero(),
            z: Fp::zero(),
        }
    }

    fn is_identity(&self) -> Choice {
        self.z.is_zero()
    }

    fn is_on_curve(&self) -> Choice {
        // Y^2 = X^3 + AX(Z^4) + b(Z^6)
        // Y^2 - (X^2 + A(Z^4))X = b(Z^6)

        let z2 = self.z.square();
        let z4 = z2.square();
        let z6 = z4 * z2;
        (self.y.square() - (self.x.square() + Ep::curve_constant_a() * z4) * self.x)
            .ct_eq(&(z6 * Ep::curve_constant_b()))
            | self.z.is_zero()
    }
}

#[cfg(test)]
mod tests {
    use super::{Ep, EpAffine};

    #[test]
    fn test_is_on_curve() {
        let affine_generator = EpAffine::generator();
        let projective_generator = Ep::generator();

        assert_eq!(affine_generator.is_on_curve().unwrap_u8(), 1);
        assert_eq!(projective_generator.is_on_curve().unwrap_u8(), 1);
    }
}
