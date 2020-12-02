# Run on CRI
# See more details in README.md
# For simplicity, we work on chr21 and chr22

module load gcc/6.2.0; module load bgen

conda activate ukb_idp

plink2_exec=/gpfs/data/im-lab/nas40t2/yanyul/softwares/plink2
outdir=/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/imagexcan_test_run
mkdir -p $outdir

# step0: subset genotype to pre-select 
# select the first 500 individuals in IDP cohort
# make a copy of plink BED genotypes.
# make a copy of BGEN genotypes.
input_geno_prefix=/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/subset_genotypes/IDP_HM3_finalPheno.chr
indiv_list=$outdir/tmp_indiv_list.txt
output_geno_prefix=$outdir/geno_for_test.chr
chr=22
cat $input_geno_prefix$chr.fam | head -n 500 > $indiv_list
for chr in `seq 21 22`
do
  $plink2_exec --bfile $input_geno_prefix$chr --keep $indiv_list --make-bed --out $output_geno_prefix$chr
  $plink2_exec --bfile $output_geno_prefix$chr --export bgen-1.2 'ref-first' --out $output_geno_prefix$chr
  bgenix -g $output_geno_prefix$chr.bgen -index -clobber
done 

# step1: simulate IDP model weights
# instead of simulating, we use the t1 weights directly.
input_step1=/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/gw_ridge/gw_ridge_beta.t1.default_theta_g_fold_5_5.parquet
output_step1=$outdir/idp_weights.parquet
if [[ ! -f $output_step1 ]]
then
  python simulate_idp_weights.py \
    $input_step1 \
    $output_step1
fi

# step2: predict IDP values
# this part needs qsub since it takes long.
input_step2=$outdir/idp_weights.parquet
output_step2=$outdir/predicted_idp.parquet
if[[ ! -f $output_step2 ]]
then
  bash pred_idp.qsub 
fi

# step3: simulate phenotype
input_step3=$outdir/predicted_idp.parquet
output_step3=$outdir/phenotype.csv
output_step3_marginal=$outdir/phenotype_marginal.yaml
output_step3_susie=$outdir/phenotype_susie.yaml
if[[ ! -f $input_step3 ]]
then
  echo "We have not run through step 2. So, cannot proceed."
  exit 1
fi
if[[ ! -f $output_step3 ]]
then
  python simulate_phenotype.py \
    $input_step3 \
    $output_step3 \
    $output_step3_marginal \
    $output_step3_susie
fi

# step4: run imagexcan
input_step4=$outdir/phenotype.csv
output_step4_marginal=$outdir/imagexcan_marginal.csv
output_step4_susie=$outdir/imagexcan_susie.csv
if[[ ! -f $output_step4 ]]
then
  python /gpfs/data/im-lab/nas40t2/yanyul/GitHub/ukb_idp_genetic_arch/methods/imagexcan/run_imagexcan.py \
    --phenotype_table $input_step4 indiv \
    --phenotype_yaml $output_step3_marginal \
    --idp_table $input_step3 indiv \
    --output $output_step4_marginal
  python /gpfs/data/im-lab/nas40t2/yanyul/GitHub/ukb_idp_genetic_arch/methods/imagexcan/run_imagexcan.py \
    --phenotype_table $input_step4 indiv \
    --phenotype_yaml $output_step3_susie \
    --idp_table $input_step3 indiv \
    --output $output_step4_susie
fi

# step5: run gwas on phenotype
output_step5=$outdir/gwas_phenotype.tsv.gz
if[[ ! -f $output_step5 ]]
then
  python /gpfs/data/im-lab/nas40t2/yanyul/GitHub/misc-tools/gw_qtl/run_trans_qtl.py \
    --geno_bed_prefix /gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/subset_genotypes/IDP_HM3_finalPheno.chr21 \
    --phenotype_table $input_step4 indiv \
    --output_prefix $outdir/gwas_phenotype.chr21 \
    --map_trans_params /gpfs/data/im-lab/nas40t2/yanyul/GitHub/misc-tools/gw_qtl/config.test.yaml
    
  python /gpfs/data/im-lab/nas40t2/yanyul/GitHub/misc-tools/gw_qtl/run_trans_qtl.py \
    --geno_bed_prefix /gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/subset_genotypes/IDP_HM3_finalPheno.chr22 \
    --phenotype_table $input_step4 indiv \
    --output_prefix $outdir/gwas_phenotype.chr22 \
    --map_trans_params /gpfs/data/im-lab/nas40t2/yanyul/GitHub/misc-tools/gw_qtl/config.test.yaml 
    
  python postprocess_gwas.py \
    $outdir/gwas_phenotype.chr21.parquet \
    $outdir/gwas_phenotype.chr22.parquet \
    $output_step5
fi

# step6: compute genotype covariance
output_step6_1=$outdir/geno_covar.chr21.naive.h5
output_step6_2=$outdir/geno_covar.chr22.naive.h5
if[[ ! -f $output_step6_1 ]]
then
  python /gpfs/data/im-lab/nas40t2/yanyul/GitHub/ukb_idp_genetic_arch/methods/simagexcan/build_genotype_covariance.py \
    --genotype_bed /gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/subset_genotypes/IDP_HM3_finalPheno.chr21.bed \
    --mode naive \
    --nbatch 50 \
    --output_prefix $outdir/geno_covar.chr21
fi
if[[ ! -f $output_step6_2 ]]
then
  python /gpfs/data/im-lab/nas40t2/yanyul/GitHub/ukb_idp_genetic_arch/methods/simagexcan/build_genotype_covariance.py \
    --genotype_bed /gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/subset_genotypes/IDP_HM3_finalPheno.chr22.bed \
    --mode naive \
    --nbatch 50 \
    --output_prefix $outdir/geno_covar.chr22
fi

# step7: s-imagexcan
output_step7=$outdir/simagexcan.csv
if[[ ! -f $output_step7 ]]
then
  python /gpfs/data/im-lab/nas40t2/yanyul/GitHub/ukb_idp_genetic_arch/methods/simagexcan/run_simagexcan.py \
    --genotype_covariance $outdir/geno_covar.chr{chr_num}.naive.h5 \
    --gwas \
    --idp_weight $output_step1 \
    --output $output_step7
fi


