Very important note:

In `pandas-plink`, it load `a0` (first column in BIM file) as reference allele, and `a1` (second column in BIM file) as alternative allele. See [here](https://pandas-plink.readthedocs.io/en/latest/usage.html#genotype) for more information.

This is **the opposite** of what is actually defined in PLINK BIM format. In other word, in PLINK BIM, usually the first column is alternative allele (dosage allele) and the second column is for the reference allele. `tensorqtl` follows this PLINK convention and flips `pandas-plink` internally (see [here](https://github.com/broadinstitute/tensorqtl/blob/master/tensorqtl/genotypeio.py#L137)).

So, to parse `tensorqtl` output, we follows the PLINK convention.
