use pairing::bn256::Fq as Fp;

#[derive(Clone, Debug)]
pub(crate) struct FfTree {
    // evaluation domain same size with polynomial
    pub(crate) domain: (Vec<Fp>, Vec<Fp>),
    // factor for performing multiplication
    pub(crate) factor: Vec<((Fp, Fp), (Fp, Fp))>,
    // inverse factor for performing multiplication
    pub(crate) inv_factor: Vec<((Fp, Fp), (Fp, Fp))>,
}

impl FfTree {
    #[cfg(test)]
    pub(crate) fn get_domain(&self) -> &(Vec<Fp>, Vec<Fp>) {
        &self.domain
    }

    pub(crate) fn get_factor(&self) -> &Vec<((Fp, Fp), (Fp, Fp))> {
        &self.factor
    }

    pub(crate) fn get_inv_factor(&self) -> &Vec<((Fp, Fp), (Fp, Fp))> {
        &self.inv_factor
    }

    pub(crate) fn last_tree(s: Vec<Fp>, s_prime: Vec<Fp>) -> Self {
        FfTree {
            domain: (s, s_prime),
            factor: vec![],
            inv_factor: vec![],
        }
    }
}
