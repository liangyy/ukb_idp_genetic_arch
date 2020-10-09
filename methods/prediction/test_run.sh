bgen=/vol/bmd/data/ukbiobank/genotypes/v3/ukb_imp_chr{chr_num}_v3.bgen
bgi=/vol/bmd/yanyul/UKB/ukb_imp_bgi/ukb_imp_chr{chr_num}_v3.bgen.bgi
sample=/vol/bmd/data/ukbiobank/genotypes/v3/ukb19526_imp_chr1_v3_s487395.sample

prs_parquet=/vol/bmd/yanyul/GitHub/ukb_idp_genetic_arch/methods/gw_ridge/test_beta.parquet
tmp_prs=test_prs_weights
tmp_parquet=$tmp_prs.parquet
tmp_txt=$tmp_prs.txt

python test_subset_parquet.py --input $prs_parquet --output_prefix $tmp_prs

nthread=4
ukb_imp_reader_path=/vol/bmd/yanyul/GitHub/misc-tools/bgen_io

output_prefix=test_prs

# conda activate ukb_idp

python run_prs.py \
  --ukb_bgen_pattern $bgen \
  --ukb_bgi_pattern $bgi \
  --ukb_sample_file $sample \
  --prs_parquet $tmp_parquet \
  --nthread $nthread \
  --output $output_prefix.parquet \
  --ukb_imp_reader_path $ukb_imp_reader_path
  

/vol/bmd/yanyul/softwares/plink2 --bfile /vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/subset_genotypes/IDP_HM3_finalPheno.merged_all --score $tmp_txt 1 3 header-read cols=scoresums,denom --score-col-nums 5 6 7 --out $output_prefix

python test_check.py
