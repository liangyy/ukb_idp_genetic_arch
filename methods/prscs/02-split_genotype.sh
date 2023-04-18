#! /usr/bin/env bash

#PBS -N split
#PBS -S /bin/bash
#PBS -l walltime=6:00:00
#PBS -l mem=4gb
#PBS -l nodes=1:ppn=1

# SPECIFY LOGGING BEHAVIOR

#PBS -o logs/split/${PBS_JOBNAME}.${PBS_JOBID}.log
#PBS -e logs/split/${PBS_JOBNAME}.${PBS_JOBID}.err

module load gcc/6.2.0
module load plink/2.0

for chr in {1..22}
do
  # train
  plink2 \
    --bfile /gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/subset_genotypes/IDP_HM3_finalPheno.chr${chr} \
    --keep /gpfs/data/im-lab/nas40t2/festus/brainxcan/train.ind.txt \
    --make-bed \
    --out /gpfs/data/im-lab/nas40t2/festus/brainxcan/genotype/train/IDP_HM3_finalPheno.chr${chr}
    
  # test
  plink2 \
    --bfile /gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/subset_genotypes/IDP_HM3_finalPheno.chr${chr} \
    --keep /gpfs/data/im-lab/nas40t2/festus/brainxcan/test.ind.txt \
    --make-bed \
    --out /gpfs/data/im-lab/nas40t2/festus/brainxcan/genotype/test/IDP_HM3_finalPheno.chr${chr}
  
done  

# plink 1.9
plink \
--merge-list mbeds.txt \
--recode vcf-iid \
--out /gpfs/data/im-lab/nas40t2/festus/brainxcan/genotype/test/IDP_HM3_finalPheno.merged_all

# plink 2.0
module load bgen/1.1.3
plink2 \
--bfile IDP_HM3_finalPheno.merged_all \
--export bgen-1.3 \
--out IDP_HM3_finalPheno.merged_all.chr1-22_v3

bgenix -index -g IDP_HM3_finalPheno.merged_all.chr1-22_v3.bgen 
