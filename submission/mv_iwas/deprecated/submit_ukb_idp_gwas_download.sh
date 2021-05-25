#!/bin/sh

outdir=/project2/haky/yanyul/data/mv_iwas/ukb_idp_gwas
runscript=/project2/haky/yanyul/GitHub/ukb_idp_genetic_arch/submission/mv_iwas/prep_ukb_idp_gwas.sh

idp_list=/project2/haky/yanyul/GitHub/ukb_idp_genetic_arch/submission/mv_iwas/ukb_idp_list.txt

position=/project2/haky/yanyul/GitHub/ukb_idp_genetic_arch/submission/mv_iwas/ukb_idp_gwas.positions.txt
if [[ ! -f $position ]]
then
  bash /project2/haky/yanyul/GitHub/ukb_idp_genetic_arch/submission/mv_iwas/prep_ukb_idp_gwas_position.sh
fi

cd $outdir 

while IFS=$'\t' read -r kk url fn idptag rem
do
    bash $runscript $url $fn $idptag $outdir $position download 
done < $idp_list


