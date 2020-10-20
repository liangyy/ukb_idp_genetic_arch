#PBS -S /bin/bash
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=32gb
#PBS -e gen_grm.err
#PBS -o gen_grm.out
#PBS -N gen_grm


cd /gpfs/data/im-lab/nas40t2/yanyul/GitHub/ukb_idp_genetic_arch/misc_data

genotype=/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/subset_genotypes/IDP_HM3_finalPheno.merged_all
plink2_exec=/gpfs/data/im-lab/nas40t2/yanyul/softwares/plink2

# output
out_prefix=/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/subset_genotypes/IDP_HM3_finalPheno.merged_all.maf_gt_0.05

$plink2_exec \
  --pfile $genotype 'vzs' \
  --make-grm-list \
  --out $out_prefix \
  --maf 0.05 \
  --threads 4 \
  --memory 30000

gzip $out_prefix.grm 

