phenotype=/vol/bmd/meliao/data/psychiatric_trait_phenotypes/2020-04-10_collected-phenotypes.txt
covariate=/vol/bmd/data/ukbiobank/psychiatric_traits/2019-12-17_psychiatric-trait_covariates.csv
idp=/vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/prediction/pred_idp.gw_ridge.parquet
indiv_list=/vol/bmd/yanyul/UKB/predicted_expression_tf2/British.txt

tmp_indiv=tmp_indiv.txt
cat $indiv_list | tail -n +2 | cut -f 1 -d ' '  > $tmp_indiv

thisdir=`pwd`

# conda activate ukb_idp
export PYTHONPATH=/vol/bmd/yanyul/GitHub/misc-tools/pyutil
python ../run_imagexcan.py \
  --covariate_table $covariate eid \
  --phenotype_table $phenotype eid \
  --covariate_yaml covar.yaml \
  --phenotype_yaml pheno.yaml \
  --individual_list $tmp_indiv \
  --idp_table $idp indiv \
  --first_30_idp \
  --output output.test_run.csv



