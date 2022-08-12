//! This curve is used for generating ecfft params
//!
//! y^2 = x^3 + x + 0x34524f71a21a7096c6ad51a6a7fe7a43d76d8e4277cb9048edc87ab655e55142
//! over past curve Fp p = 0x40000000000000000000000000000000224698fc094cf91b992d30ed00000001
//! a = 1
//! b = 23665697887148517506426806798051226694671519983424102823343279587811911881026
//! n = 2^12 | #E(Fp)
//!
//! G ⊂ E(Fp)
//! #G = n
//!
//! These params allow us to evaluate `n` degree polynomials

use core::fmt::{self, Debug};
use core::ops::{Add, Mul};

use super::behave::{
    curve_affine_coordinate_method, curve_constant_params, curve_projective_arithmetic,
    curve_projective_coordinate_method,
};
use pairing::bn256::Fr as Fp;
use pairing::group::ff::{Field, PrimeField};
use subtle::{Choice, ConditionallySelectable, ConstantTimeEq};

const CURVE_A: Fp = Fp::one();

const CURVE_B: Fp = Fp::from_raw([
    0x5dcdee14b5ed61a0,
    0x35df7da06ba32982,
    0x3bfb29b83daa1fd1,
    0xc6871bc29d46163,
]);

const GENERATOR_X: Fp = Fp::from_raw([
    0xa6cbe58e064a6b0c,
    0xe067e9f0c64a3125,
    0xb2875bc26ec52bee,
    0xcff00b8abd9119b,
]);

const GENERATOR_Y: Fp = Fp::from_raw([
    0xc0b0b0de539981b4,
    0x2d43482567de0148,
    0x0ac0f5b3615a7946,
    0x2ef5a31d89db82f6,
]);

const REPRESENTATIVE_X: Fp = Fp::from_raw([
    0x3013ff91324dd873,
    0x497c46cd5f370537,
    0x3390169392010587,
    0x19d9b63388a8b653,
]);

const REPRESENTATIVE_Y: Fp = Fp::from_raw([
    0x9d3ee38d780e552e,
    0x867ccec0c36810ab,
    0xc9d8e37e568729c1,
    0x737e1ab7c2b3a75,
]);

curve_affine_coordinate_method!(EpAffine, Fp, CURVE_A, CURVE_B);
curve_projective_coordinate_method!(Ep, Fp, CURVE_A, CURVE_B);

curve_projective_arithmetic!(Ep, EpAffine, Fp);

curve_constant_params!(
    Ep,
    EpAffine,
    Fp,
    GENERATOR_X,
    GENERATOR_Y,
    REPRESENTATIVE_X,
    REPRESENTATIVE_Y
);

#[cfg(test)]
mod tests {
    use super::{Ep, EpAffine, Fp};
    use pairing::group::ff::Field;
    use rand_core::OsRng;

    #[test]
    fn test_const_points_is_on_curve() {
        let affine_generator = EpAffine::generator();
        let projective_generator = Ep::generator();

        let affine_representative = EpAffine::representative();
        let projective_representative = Ep::representative();

        assert_eq!(affine_generator.is_on_curve().unwrap_u8(), 1);
        assert_eq!(projective_generator.is_on_curve().unwrap_u8(), 1);

        assert_eq!(affine_representative.is_on_curve().unwrap_u8(), 1);
        assert_eq!(projective_representative.is_on_curve().unwrap_u8(), 1);
    }

    #[test]
    fn test_add_points() {
        let projective_representative = Ep::representative();
        let projective_generator = Ep::generator();

        let add = projective_generator + projective_representative;
        assert_eq!(add.is_on_curve().unwrap_u8(), 1);
    }

    #[test]
    fn test_double_point() {
        let projective_representative = Ep::representative();

        let double = projective_representative.double();
        assert_eq!(double.is_on_curve().unwrap_u8(), 1);
    }

    #[test]
    fn test_subgroup_order() {
        let projective_generator = Ep::generator();
        let order = 1 << 14;
        let identity = projective_generator * Fp::from_raw([order, 0, 0, 0]);

        assert_eq!(identity.is_identity().unwrap_u8(), 1);
    }

    #[test]
    fn test_to_affine() {
        let projective_generator = Ep::generator();
        let random_point = projective_generator * Fp::random(OsRng);
        let affine_point = random_point.to_affine();

        assert_eq!(affine_point.is_on_curve().unwrap_u8(), 1);
    }
}
