# args1: idp matrix
# args2: outdir
# args3: nametag

# environment setup
source /vol/bmd/yanyul/miniconda3/etc/profile.d/conda.sh
conda activate ukb_idp

cd /vol/bmd/yanyul/GitHub/ukb_idp_genetic_arch/methods/heritability

runscript=run_pyemma.sh
grm=/vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/subset_genotypes/IDP_HM3_finalPheno.merged_all.maf_gt_0.05
idp=$1  # /vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/regress_out_idp_pcs/2020-05-18_final-phenotypes.cleaned_up_T1.parquet
outdir=$2  # /vol/bmd/yanyul/UKB/ukb_idp_genetic_arch/heritability_2nd_round

nametag=$3  # pyemma.idp_t1


screen -dmS $nametag bash -c "bash $runscript $idp $grm $outdir/$nametag"