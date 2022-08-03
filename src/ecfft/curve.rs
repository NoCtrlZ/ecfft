//! This curve is used for generating ecfft params
//!
//! y^2 = x^3 + x + 0x34524f71a21a7096c6ad51a6a7fe7a43d76d8e4277cb9048edc87ab655e55142
//! over past curve Fp p = 0x40000000000000000000000000000000224698fc094cf91b992d30ed00000001
//! a = 1
//! b = 23665697887148517506426806798051226694671519983424102823343279587811911881026
//! n = 2^12 | #E(Fp)

use core::fmt::{self, Debug};
use core::ops::{Add, Mul};

use super::behave::{
    curve_affine_coordinate_method, curve_projective_arithmetic, curve_projective_coordinate_method,
};
use ff::{Field, PrimeField};
use pasta_curves::Fp;
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};

curve_affine_coordinate_method!(
    EpAffine,
    Fp,
    [1, 0, 0, 0],
    [
        0xedc87ab655e55142,
        0xd76d8e4277cb9048,
        0xc6ad51a6a7fe7a43,
        0x34524f71a21a7096,
    ]
);

impl EpAffine {
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

    fn subgroup_generator() -> Self {
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
        }
    }

    fn representative() -> Self {
        Self {
            x: Fp::from_raw([
                0x5b1575e3d64364c7,
                0x8f4b66263e159043,
                0xb8787be360f5270b,
                0x1ac2d405d53f2333,
            ]),
            y: Fp::from_raw([
                0x641bc6c587cf1cb9,
                0x6082deb9cca88fef,
                0x802db2846bbce985,
                0x3488bbb09671330,
            ]),
        }
    }

    fn to_curve(&self) -> Ep {
        Ep {
            x: self.x,
            y: self.y,
            z: Fp::conditional_select(&Fp::one(), &Fp::zero(), self.is_identity()),
        }
    }
}

curve_projective_coordinate_method!(
    Ep,
    Fp,
    [1, 0, 0, 0],
    [
        0xedc87ab655e55142,
        0xd76d8e4277cb9048,
        0xc6ad51a6a7fe7a43,
        0x34524f71a21a7096,
    ]
);
curve_projective_arithmetic!(Ep, EpAffine, Fp);

impl Ep {
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

    fn subgroup_generator() -> Self {
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

    fn representative() -> Self {
        Self {
            x: Fp::from_raw([
                0x5b1575e3d64364c7,
                0x8f4b66263e159043,
                0xb8787be360f5270b,
                0x1ac2d405d53f2333,
            ]),
            y: Fp::from_raw([
                0x641bc6c587cf1cb9,
                0x6082deb9cca88fef,
                0x802db2846bbce985,
                0x3488bbb09671330,
            ]),
            z: Fp::from_raw([1, 0, 0, 0]),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{Ep, EpAffine};

    #[test]
    fn test_const_points_is_on_curve() {
        let affine_generator = EpAffine::generator();
        let projective_generator = Ep::generator();

        let affine_subgroup_generator = EpAffine::subgroup_generator();
        let projective_subgroup_generator = Ep::subgroup_generator();

        let affine_representative = EpAffine::representative();
        let projective_representative = Ep::representative();

        assert_eq!(affine_generator.is_on_curve().unwrap_u8(), 1);
        assert_eq!(projective_generator.is_on_curve().unwrap_u8(), 1);

        assert_eq!(affine_subgroup_generator.is_on_curve().unwrap_u8(), 1);
        assert_eq!(projective_subgroup_generator.is_on_curve().unwrap_u8(), 1);

        assert_eq!(affine_representative.is_on_curve().unwrap_u8(), 1);
        assert_eq!(projective_representative.is_on_curve().unwrap_u8(), 1);
    }

    #[test]
    fn test_add_points() {
        let projective_representative = Ep::representative();
        let projective_subgroup_generator = Ep::subgroup_generator();

        let add = projective_subgroup_generator + projective_representative;
        assert_eq!(add.is_on_curve().unwrap_u8(), 1);
    }

    #[test]
    fn test_double_point() {
        let projective_representative = Ep::representative();

        let double = projective_representative.double();
        assert_eq!(double.is_on_curve().unwrap_u8(), 1);
    }
}
