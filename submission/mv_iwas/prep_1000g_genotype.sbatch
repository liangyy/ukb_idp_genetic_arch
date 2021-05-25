#!/bin/sh

#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --mem-per-cpu=8G
#SBATCH --nodes=2

module load plink
module load htslib
module load bcftools

outdir=/project2/haky/yanyul/data/mv_iwas/1000g
runscript=/project2/haky/yanyul/GitHub/ukb_idp_genetic_arch/submission/mv_iwas/prep_1000g_genotype.sh
eurlist=/project2/haky/yanyul/GitHub/ukb_idp_genetic_arch/submission/mv_iwas/EUR.txt

cd $outdir 

bash $runscript $CHR $outdir $eurlist 

