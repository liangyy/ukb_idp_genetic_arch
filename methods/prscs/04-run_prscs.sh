#! /usr/bin/env bash

# t1 region
###path.t=/gpfs/data/im-lab/nas40t2/yanyul/data/brainxcan_data/idp_gwas/residual.t1.chr1 #original
###path.t=/scratch/fnyasimi1/ukb_idp/idp_gwas_4th/trans_qtl.fourth_round.t1_w_pc.chr1
for gg in $(ls /scratch/fnyasimi1/ukb_idp/idp_gwas_4th/trans_qtl.fourth_round.t1_w_pc.chr1/*.parquet)
do
    pre=$(basename ${gg} .parquet)
    [[ ! -f results.all/t1/${pre}.tsv ]] && echo "Submitting ${pre}" && qsub -v region=t1,gwas_ss=${gg} prscs_runner.pbs 
done


# dmri region
###path.d=/gpfs/data/im-lab/nas40t2/yanyul/data/brainxcan_data/idp_gwas/residual.t1.chr1 #original
###path.d=/scratch/fnyasimi1/ukb_idp/idp_gwas_4th/trans_qtl.fourth_round.dmri_w_pc.chr1
for gg in $(ls /scratch/fnyasimi1/ukb_idp/idp_gwas_4th/trans_qtl.fourth_round.dmri_w_pc.chr1/*.parquet)
do
    pre=$(basename ${gg} .parquet)
    [[ ! -f results.all/dmri/${pre}.tsv ]] && echo "Submitting ${pre}" && qsub -v region=dmri,gwas_ss=${gg} prscs_runner.pbs 
done


