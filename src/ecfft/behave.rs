macro_rules! curve_projective_coordinate_method {
    ($curve:ident, $field:ident, $a:expr, $b:expr) => {
        #[derive(Copy, Clone, Debug)]
        pub(crate) struct $curve {
            x: $field,
            y: $field,
            z: $field,
        }

        impl $curve {
            const fn curve_constant_a() -> $field {
                $field::from_raw($a)
            }

            const fn curve_constant_b() -> $field {
                $field::from_raw($b)
            }

            fn is_identity(&self) -> Choice {
                self.z.is_zero()
            }

            fn identity() -> $curve {
                $curve {
                    x: $field::zero(),
                    y: $field::zero(),
                    z: $field::zero(),
                }
            }

            fn is_on_curve(&self) -> Choice {
                let z2 = self.z.square();
                let z4 = z2.square();
                let z6 = z4 * z2;
                (self.y.square() - (self.x.square() + Ep::curve_constant_a() * z4) * self.x)
                    .ct_eq(&(z6 * Ep::curve_constant_b()))
                    | self.z.is_zero()
            }
        }

        impl PartialEq for $curve {
            fn eq(&self, other: &Self) -> bool {
                self.ct_eq(other).into()
            }
        }

        impl Eq for $curve {}

        impl ConstantTimeEq for $curve {
            fn ct_eq(&self, other: &Self) -> Choice {
                // Is (xz^2, yz^3, z) equal to (x'z'^2, yz'^3, z') when converted to affine?

                let z = other.z.square();
                let x1 = self.x * z;
                let z = z * other.z;
                let y1 = self.y * z;
                let z = self.z.square();
                let x2 = other.x * z;
                let z = z * self.z;
                let y2 = other.y * z;

                let self_is_zero = self.is_identity();
                let other_is_zero = other.is_identity();

                (self_is_zero & other_is_zero) // Both point at infinity
                            | ((!self_is_zero) & (!other_is_zero) & x1.ct_eq(&x2) & y1.ct_eq(&y2))
                // Neither point at infinity, coordinates are the same
            }
        }

        impl ConditionallySelectable for $curve {
            fn conditional_select(a: &Self, b: &Self, choice: Choice) -> Self {
                $curve {
                    x: $field::conditional_select(&a.x, &b.x, choice),
                    y: $field::conditional_select(&a.y, &b.y, choice),
                    z: $field::conditional_select(&a.z, &b.z, choice),
                }
            }
        }
    };
}

macro_rules! curve_projective_arithmetic {
    ($curve:ident, $curve_affine:ident, $field:ident) => {
        impl $curve {
            fn to_affine(&self) -> $curve_affine {
                let zinv = self.z.invert().unwrap_or($field::zero());
                let zinv2 = zinv.square();
                let x = self.x * zinv2;
                let zinv3 = zinv2 * zinv;
                let y = self.y * zinv3;

                $curve_affine { x, y }
            }

            fn double(&self) -> Self {
                // https://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#doubling-dbl-2007-bl

                let xx = self.x.square();
                let yy = self.y.square();
                let yyyy = yy.square();
                let zz = self.z.square();
                let a = self.x + yy;
                let b = a.square() - xx - yyyy;
                let s = b + b;
                let c = xx + xx + xx;
                let d = zz.square() * Self::curve_constant_a();
                let m = c + d;
                let e = s + s;
                let t = m.square() - e;
                let x3 = t;
                let f = s - t;
                let l = yyyy.double().double().double();
                let y3 = m * f - l;
                let n = self.y + self.z;
                let z3 = n.square() - yy - zz;

                let tmp = $curve {
                    x: x3,
                    y: y3,
                    z: z3,
                };

                $curve::conditional_select(&tmp, &$curve::identity(), self.is_identity())
            }
        }

        impl<'a, 'b> Add<&'a Ep> for &'b $curve {
            type Output = $curve;

            fn add(self, rhs: &'a $curve) -> $curve {
                if bool::from(self.is_identity()) {
                    *rhs
                } else if bool::from(rhs.is_identity()) {
                    *self
                } else {
                    let z1z1 = self.z.square();
                    let z2z2 = rhs.z.square();
                    let u1 = self.x * z2z2;
                    let u2 = rhs.x * z1z1;
                    let s1 = self.y * z2z2 * rhs.z;
                    let s2 = rhs.y * z1z1 * self.z;

                    if u1 == u2 {
                        if s1 == s2 {
                            self.double()
                        } else {
                            $curve::identity()
                        }
                    } else {
                        let h = u2 - u1;
                        let i = (h + h).square();
                        let j = h * i;
                        let r = s2 - s1;
                        let r = r + r;
                        let v = u1 * i;
                        let x3 = r.square() - j - v - v;
                        let s1 = s1 * j;
                        let s1 = s1 + s1;
                        let y3 = r * (v - x3) - s1;
                        let z3 = (self.z + rhs.z).square() - z1z1 - z2z2;
                        let z3 = z3 * h;

                        $curve {
                            x: x3,
                            y: y3,
                            z: z3,
                        }
                    }
                }
            }
        }

        impl<'a, 'b> Add<&'a $curve_affine> for &'b $curve {
            type Output = $curve;

            fn add(self, rhs: &'a $curve_affine) -> $curve {
                if bool::from(self.is_identity()) {
                    rhs.to_curve()
                } else if bool::from(rhs.is_identity()) {
                    *self
                } else {
                    let z1z1 = self.z.square();
                    let u2 = rhs.x * z1z1;
                    let s2 = rhs.y * z1z1 * self.z;

                    if self.x == u2 {
                        if self.y == s2 {
                            self.double()
                        } else {
                            $curve::identity()
                        }
                    } else {
                        let h = u2 - self.x;
                        let hh = h.square();
                        let i = hh + hh;
                        let i = i + i;
                        let j = h * i;
                        let r = s2 - self.y;
                        let r = r + r;
                        let v = self.x * i;
                        let x3 = r.square() - j - v - v;
                        let j = self.y * j;
                        let j = j + j;
                        let y3 = r * (v - x3) - j;
                        let z3 = (self.z + h).square() - z1z1 - hh;

                        $curve {
                            x: x3,
                            y: y3,
                            z: z3,
                        }
                    }
                }
            }
        }

        impl<'a, 'b> Add<&'a $curve> for &'b $curve_affine {
            type Output = $curve;

            fn add(self, rhs: &'a $curve) -> $curve {
                rhs + self
            }
        }

        impl<'a, 'b> Add<&'a $curve_affine> for &'b $curve_affine {
            type Output = $curve;

            fn add(self, rhs: &'a $curve_affine) -> $curve {
                if bool::from(self.is_identity()) {
                    rhs.to_curve()
                } else if bool::from(rhs.is_identity()) {
                    self.to_curve()
                } else {
                    if self.x == rhs.x {
                        if self.y == rhs.y {
                            self.to_curve().double()
                        } else {
                            $curve::identity()
                        }
                    } else {
                        let h = rhs.x - self.x;
                        let hh = h.square();
                        let i = hh + hh;
                        let i = i + i;
                        let j = h * i;
                        let r = rhs.y - self.y;
                        let r = r + r;
                        let v = self.x * i;
                        let x3 = r.square() - j - v - v;
                        let j = self.y * j;
                        let j = j + j;
                        let y3 = r * (v - x3) - j;
                        let z3 = h + h;

                        $curve {
                            x: x3,
                            y: y3,
                            z: z3,
                        }
                    }
                }
            }
        }

        impl<'b> Add<&'b $curve> for $curve {
            type Output = $curve;

            #[inline]
            fn add(self, rhs: &'b $curve) -> $curve {
                &self + rhs
            }
        }

        impl<'a> Add<$curve> for &'a $curve {
            type Output = $curve;

            #[inline]
            fn add(self, rhs: $curve) -> $curve {
                self + &rhs
            }
        }

        impl Add<$curve> for $curve {
            type Output = $curve;

            #[inline]
            fn add(self, rhs: $curve) -> $curve {
                &self + &rhs
            }
        }

        impl<'a, 'b> Mul<&'b $field> for &'a $curve {
            type Output = $curve;

            fn mul(self, other: &'b $field) -> Self::Output {
                let mut acc = $curve::identity();

                for bit in other
                    .to_repr()
                    .iter()
                    .rev()
                    .flat_map(|byte| (0..8).rev().map(move |i| Choice::from((byte >> i) & 1u8)))
                    .skip(1)
                {
                    acc = acc.double();
                    acc = $curve::conditional_select(&acc, &(acc + self), bit);
                }

                acc
            }
        }

        impl<'b> Mul<&'b $field> for $curve {
            type Output = $curve;

            #[inline]
            fn mul(self, rhs: &'b $field) -> $curve {
                &self * rhs
            }
        }

        impl<'a> Mul<$field> for &'a $curve {
            type Output = $curve;

            #[inline]
            fn mul(self, rhs: $field) -> $curve {
                self * &rhs
            }
        }

        impl Mul<$field> for $curve {
            type Output = $curve;

            #[inline]
            fn mul(self, rhs: $field) -> $curve {
                &self * &rhs
            }
        }
    };
}

macro_rules! curve_affine_coordinate_method {
    ($curve_affine:ident, $field:ident, $a:expr, $b:expr) => {
        #[derive(Clone)]
        pub(crate) struct $curve_affine {
            x: $field,
            y: $field,
        }

        impl $curve_affine {
            const fn curve_constant_a() -> $field {
                $field::from_raw($a)
            }

            const fn curve_constant_b() -> $field {
                $field::from_raw($b)
            }

            fn is_on_curve(&self) -> Choice {
                (self.y.square() - (self.x.square() + $curve_affine::curve_constant_a()) * self.x)
                    .ct_eq(&$curve_affine::curve_constant_b())
                    | self.is_identity()
            }

            fn is_identity(&self) -> Choice {
                self.x.is_zero() & self.y.is_zero()
            }
        }

        impl fmt::Debug for $curve_affine {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
                if self.is_identity().into() {
                    write!(f, "Infinity")
                } else {
                    write!(f, "({:?}, {:?})", self.x, self.y)
                }
            }
        }
    };
}

macro_rules! curve_constant_params {
    ($curve:ident, $curve_affine:ident, $field:ident, $generator_x:ident, $generator_y:ident, $subgroup_generator_x:ident, $subgroup_generator_y:ident, $representative_x:ident, $representative_y:ident) => {
        impl $curve {
            fn generator() -> Self {
                Self {
                    x: $generator_x,
                    y: $generator_y,
                    z: Fp::from_raw([1, 0, 0, 0]),
                }
            }

            fn subgroup_generator() -> Self {
                Self {
                    x: $subgroup_generator_x,
                    y: $subgroup_generator_y,
                    z: Fp::from_raw([1, 0, 0, 0]),
                }
            }

            fn representative() -> Self {
                Self {
                    x: $representative_x,
                    y: $representative_y,
                    z: Fp::from_raw([1, 0, 0, 0]),
                }
            }
        }

        impl $curve_affine {
            fn generator() -> Self {
                $curve_affine {
                    x: $generator_x,
                    y: $generator_y,
                }
            }

            fn subgroup_generator() -> Self {
                $curve_affine {
                    x: $subgroup_generator_x,
                    y: $subgroup_generator_y,
                }
            }

            fn representative() -> Self {
                $curve_affine {
                    x: $representative_x,
                    y: $representative_y,
                }
            }

            fn to_curve(&self) -> $curve {
                $curve {
                    x: self.x,
                    y: self.y,
                    z: $field::conditional_select(
                        &$field::one(),
                        &$field::zero(),
                        self.is_identity(),
                    ),
                }
            }
        }
    };
}

pub(crate) use {
    curve_affine_coordinate_method, curve_constant_params, curve_projective_arithmetic,
    curve_projective_coordinate_method,
};
