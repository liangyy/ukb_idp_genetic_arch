
Here I implement the ridge regression for genome-wide SNPs.

**Dependencies**: Tested on `pandas=1.1.1`, `numpy=1.19.1`, and `fastparquet=0.4.1`.

The script is designed to use the genotype in PLINK binary PED format.
And it outputs the cross-validated prediction performance of the ridge regression, where the hyper-parameter is determined by a nested round of cross-validation within training folds.

Here we follow the BLUP based formula and we assume that the GRM can be fit into memory. 
If this is not the case, consider reduce the precision of the floating number in GRM. 
For instance, reducing from `float64` to `float32` will reduce the memory need by a half. 

The script to test the performance on a grid of `theta_g` where we assign `theta_g` weight to GRM and `1 - theta_g` weight to the identity component. 
`theta_g` varies between 0 and 1 by construction.
The relation between `theta_g` and `sigma2` (in vanilla ridge regression formula) is `sigma2 = (1 - theta_g) M / theta_g` where `M` is the number of SNPs.  

Testing command:

```
python run_gw_ridge.py \
  --geno_bed_pattern /vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/subset_genotypes/IDP_HM3_finalPheno.chr{chr_num}.bed \
  --phenotype_parquet /vol/bmd/meliao/data/idp_phenotypes/2020-05-18_final-phenotypes.parquet \
  --nfold 5 5 \
  --first_n_indiv 500 \
  --output test.tsv.gz
```

```
python run_gw_ridge.py  \
  --geno_bed_pattern /vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/subset_genotypes/IDP_HM3_finalPheno.chr{chr_num}.bed  \
  --phenotype_parquet /vol/bmd/meliao/data/idp_phenotypes/2020-05-18_final-phenotypes.parquet \
  --first_n_indiv 500 \
  --output test_beta.parquet \
  --nfold 5 5 \
  --train_full_model 
```

Use GCTA GRM.

```
python run_gw_ridge.py \
  --gcta_grm_prefix first_500 \
  --phenotype_parquet /vol/bmd/meliao/data/idp_phenotypes/2020-05-18_final-phenotypes.parquet \
  --nfold 5 5 \
  --output test_gcta_grm.tsv.gz 
```

```
python run_gw_ridge.py \
  --gcta_grm_prefix first_500 \
  --phenotype_parquet /vol/bmd/meliao/data/idp_phenotypes/2020-05-18_final-phenotypes.parquet \
  --output test_gcta_grm_beta.parquet \
  --train_full_model \
  --nfold 5 5 \
  --geno_bed_pattern /vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/subset_genotypes/IDP_HM3_finalPheno.chr{chr_num}.bed 
```

Load one chromosome at a time.

```
python run_gw_ridge.py \
  --geno_bed_pattern first_500.bed  \
  --phenotype_parquet /vol/bmd/meliao/data/idp_phenotypes/2020-05-18_final-phenotypes.parquet \
  --nfold 5 5 \
  --output test_load_chr.tsv.gz 
```

```
python run_gw_ridge.py \
  --geno_bed_pattern first_500.bed \
  --phenotype_parquet /vol/bmd/meliao/data/idp_phenotypes/2020-05-18_final-phenotypes.parquet \
  --output test_load_chr_beta.parquet \
  --train_full_model \
  --nfold 5 5
```

Misc scripts for testing.

```
# generate prs weights for plink
import pandas as pd
e = pd.read_parquet('test_beta.parquet')
e[['snpid', 'a0', 'a1', 'IDP-25303', 'IDP-25311']].to_csv('test_beta_idp_25311_25304.tsv', sep='\t', index=False)


# compare PRS with observed IDP
import pandas as pd
f = pd.read_parquet('/vol/bmd/meliao/data/idp_phenotypes/2020-05-18_final-phenotypes.parquet')
e = pd.read_csv('test_prs.sscore', sep='\t')

import gzip, pickle
with gzip.open('test_beta.parquet.grm_cache.pkl.gz', 'rb') as ff:
    oo = pickle.load(ff)

e['#IID'] = e['#IID'].astype(str)
e = pd.merge(e, f, left_on='#IID', right_on='individual')
e = e[ e.individual.isin(oo['grm_indiv_info']) ]
import numpy as np
print(
np.corrcoef(e['IDP-25303_SUM'], e['IDP-25303']), 
np.corrcoef(e['IDP-25311_SUM'], e['IDP-25311']),
np.corrcoef(e['IDP-25303_SUM'], e['IDP-25311']), 
np.corrcoef(e['IDP-25311_SUM'], e['IDP-25303'])
)
```

