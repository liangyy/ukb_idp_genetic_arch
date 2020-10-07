Here I implement the ridge regression for genome-wide SNPs.

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
  --output test.tsv.gz \
  --snplist_to_exclude /vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/subset_genotypes/IDP_HM3_finalPheno.merged_all-merge.missnp
```

```
python run_gw_ridge.py  \
  --geno_bed_pattern /vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/subset_genotypes/IDP_HM3_finalPheno.chr{chr_num}.bed  \
  --phenotype_parquet /vol/bmd/meliao/data/idp_phenotypes/2020-05-18_final-phenotypes.parquet \
  --nfold 5 5 \
  --first_n_indiv 500 \
  --output test_beta.parquet \
  --train_full_model \
  --snplist_to_exclude /vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/subset_genotypes/IDP_HM3_finalPheno.merged_all-merge.missnp
```
