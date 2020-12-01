condadir=`which R | sed 's#bin/R##g'`
genofile=$condadir/lib/R/library/snpnet/extdata/sample
phenofile=test_files/sample.phe.tsv
phenolist=test_files/test_pheno_list.txt
nfold=5
innernfold=5
indiv_col=IID
snpnet_config=test_files/snpnet.yaml
output_prefix=test_files/output_test_training

# conda activate snpnet

# cd /vol/bmd/yanyul/GitHub/ukb_idp_genetic_arch/methods/gw_lasso

Rscript run_gw_lasso.R \
  --genotype $genofile \
  --phenotype_table $phenofile \
  --nfold $nfold \
  --inner_nfold $innernfold \
  --indiv_col $indiv_col \
  --pheno_list $phenolist \
  --snpnet_config $snpnet_config \
  --output_prefix $output_prefix \
  --mode model_training \
  > test_files/test_run_training.log 2>&1
  
