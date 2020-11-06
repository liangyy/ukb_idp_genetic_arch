runscript=../../methods/heritability/run_pyemma.sh
grm=/vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/subset_genotypes/IDP_HM3_finalPheno.merged_all.maf_gt_0.05
# example phenotype: /vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/regress_out_idp_pcs/2020-05-18_final-phenotypes.regress_out_10PCs.parquet
idp_prefix=/vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/regress_out_idp_pcs/2020-05-18_final-phenotypes.regress_out_
idp_suffix=PCs.parquet
outdir=/vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/heritability

nametag=pyemma.idp_pc
npcs='0 1 2 5 10 20 50'

for i in $npcs
do
  
  nametag_i=$nametag$i
  screen -dmS $nametag_i bash -c "bash $runscript $idp_prefix$i$idp_suffix $grm $outdir/$nametag_i"

done
