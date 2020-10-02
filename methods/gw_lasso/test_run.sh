genofile=/vol/bmd/yanyul/miniconda3/envs/snpnet/lib/R/library/snpnet/extdata/sample
phenofile=/vol/bmd/yanyul/miniconda3/envs/snpnet/lib/R/library/snpnet/extdata/sample.phe
phenolist=test_files/test_pheno_list.txt
nfold=5
innernfold=5
indiv_col=IID
snpnet_config=test_files/snpnet.yaml
output_prefix=test_files/output

# conda activate snpnet

cd /vol/bmd/yanyul/GitHub/ukb_idp_genetic_arch/methods/gw_lasso

Rscript run_gw_lasso.R \
  --genotype $genofile \
  --phenotype_table $phenofile \
  --nfold $nfold \
  --inner_nfold $innernfold \
  --indiv_col $indiv_col \
  --pheno_list $phenolist \
  --snpnet_config $snpnet_config \
  --output_prefix $output_prefix \
  > test_run.log 2>&1
  