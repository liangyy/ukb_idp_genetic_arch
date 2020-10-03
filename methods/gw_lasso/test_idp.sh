genofile=/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/test_genotypes/first1000
phenofile=/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/test_genotypes/test_two_idps.phe.txt
phenolist=/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/test_genotypes/test_pheno_list.txt
nfold=5
innernfold=5
indiv_col=individual
snpnet_config=test_files/snpnet.yaml
output_prefix=test_files/test_idp

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
  > test_idp.log 2>&1
  
