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
  
# check whether the weights generated is as expected.
# run plink2 --score and compare with the in sample prediction output by snpnet
plink2_exe=/gpfs/data/im-lab/nas40t2/yanyul/softwares/plink2
bedgeno=/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/subset_genotypes/IDP_HM3_finalPheno.merged_all
tmp_txt=QPHE.prs.txt
zcat $output_prefix.weights.tsv.gz | grep QPHE > $tmp_txt

$plink2_exe --pfile $bedgeno 'vzs' --score $tmp_txt 1 4 2 cols=scoresums,denom --out $output_prefix