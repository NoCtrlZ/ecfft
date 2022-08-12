//! This curve is used for generating ecfft params
//!
//! y^2 = x^3 + x + 5612291247948481584627780310922020304781354847659642188369727566000581075360
//! over past curve Fp p = 0x30644E72E131A029B85045B68181585D97816A916871CA8D3C208C16D87CFD47
//! a = 1
//! b = 5612291247948481584627780310922020304781354847659642188369727566000581075360
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
use pairing::bn256::Fq as Fp;
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
    0xf3c0656f8ff34e2c,
    0x75ce795513e0bf38,
    0x260d26c12f05dd81,
    0x1c30bfe9bd1623f6,
]);

const GENERATOR_Y: Fp = Fp::from_raw([
    0x1ee3da4379313853,
    0xd158cae28c7b0238,
    0xe1c47936c576ec08,
    0x2851ed59c90496af,
]);

const REPRESENTATIVE_X: Fp = Fp::from_raw([
    0x43cc2220bc2c43d9,
    0x2b4975b0074a4c19,
    0x49e9941334ba73f5,
    0x13305a0c7fac8d4e,
]);

const REPRESENTATIVE_Y: Fp = Fp::from_raw([
    0xe85005efc204f710,
    0x70a01a702609f044,
    0x58d344b7260f9ab4,
    0x1d252a06a64cfbef,
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
