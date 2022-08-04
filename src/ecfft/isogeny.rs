use ff::Field;

// isogeny structure
#[derive(Clone, Debug)]
pub struct Isogeny<F: Field> {
    pub numerator: [F; 3],
    pub denominator: [F; 2],
}

impl<F: Field> Isogeny<F> {
    pub(crate) fn evaluate(&self, x: F) -> F {
        let Isogeny {
            numerator: [a0, a1, a2],
            denominator: [b0, b1],
        } = self;
        let numerator = *a0 + *a1 * x + *a2 * x * x;
        let denominator = *b0 + *b1 * x;
        numerator * denominator.invert().unwrap()
    }
}
