runscript=../../methods/heritability/run_pyemma.sh
grm=/vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/subset_genotypes/IDP_HM3_finalPheno.merged_all.maf_gt_0.05
idp=/vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/regress_out_idp_pcs/2020-05-18_final-phenotypes.cleaned_up_T1.parquet
outdir=/vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/heritability

nametag=pyemma.idp_t1


screen -dmS $nametag bash -c "bash $runscript $idp $grm $outdir/$nametag"
