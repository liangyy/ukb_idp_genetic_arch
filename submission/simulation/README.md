We perform a simulation study to compare the performance of BrainXcan and genetic correlation.
We assume a mediation model (see [here](https://github.com/liangyy/ukb_idp_genetic_arch/blob/master/rmd/simulation_study.Rmd) for a toy example).

In this directory, we perform similar simulation at genome-wide scale.
The procedures are listed below:

1. **Split UKB IDP genotypes** into two groups: group A and B
2. **Simulate mediators and phenotypes**. For both groups,
    - Simulate mediators with dense `B` under different h2 values: `M = X B + E`
    - Simulate phenotypes with sparse `beta` under different PVE values: `Y = M beta + e`
    - For null, we simulate phenotypes with `Y = X b + e` with dense `b`
3. **Train ridge models for mediators** (group 1)
4. **Run GWAS for phenotypes and mediators** (group 1 and 2)
5. **Calculate LD scores** using the current genotype files as the reference panel (limited SNP set)
6. **Run BrainXcan and genetic correlation**
