module load gcc/6.2.0
module load gcta

genotype=/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/subset_genotypes/IDP_HM3_finalPheno.merged_all
gcta_exec=gcta
plink2_exec=/gpfs/data/im-lab/nas40t2/yanyul/softwares/plink2

# output
out_prefix=/gpfs/data/im-lab/nas40t2/yanyul/ukb_idp/subset_genotypes/IDP_HM3_finalPheno.merged_all.maf_gt_0.05

$plink2_exec \
  --pfile $genotype \
  --make-grm \
  --out $out_prefix \
  --maf 0.05 

gzip $out_prefix.grm 
