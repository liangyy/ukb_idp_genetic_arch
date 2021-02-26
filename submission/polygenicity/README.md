# Pre-compute LD matrix for S-BayesS

```
for i in `seq 1 5`
do
  qsub -v CHRNUM=$i -l mem=64gb -l walltime=72:00:00 run_sbayess_precomp_ld.qsub
done
for i in `seq 6 16`
do
  qsub -v CHRNUM=$i -l mem=32gb -l walltime=72:00:00 run_sbayess_precomp_ld.qsub
done
for i in `seq 17 22`
do
  qsub -v CHRNUM=$i -l mem=16gb -l walltime=48:00:00 run_sbayess_precomp_ld.qsub
done
``` 

# Run S-BayesS

```
dmri=/gpfs/data/im-lab/nas40t2/yanyul/GitHub/ukb_idp_genetic_arch/submission/gw_lasso/dmri_2nd_round/list_misc/dmri.original.all_covar.w_pc/pheno_list.txt
t1=/gpfs/data/im-lab/nas40t2/yanyul/GitHub/ukb_idp_genetic_arch/submission/gw_lasso/t1_2nd_round/list_misc/t1.scaled.all_covar.w_pc/pheno_list.txt

bash submit_sbayess_run.sh $dmri
bash submit_sbayess_run.sh $t1
```