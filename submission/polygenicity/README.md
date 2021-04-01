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

bash submit_sbayess_run.sh $dmri sbayess_dmri
bash submit_sbayess_run.sh $t1 sbayess_t1
```

# Submit SLD4M jobs

On midway2.

## GWAS external

```
NAME_LIST=/project2/haky/yanyul/GitHub/ukb_idp_genetic_arch/misc_data/ld4m/ld4m_external_traits.txt
MYJOB_NAME=sld4m_external
sbatch --export=NAME_LIST=$NAME_LIST,MYJOB_NAME=$MYJOB_NAME run_sld4m.sbatch
```

## GWAS independent

```
NAME_LIST=/project2/haky/yanyul/GitHub/ukb_idp_genetic_arch/misc_data/ld4m/ld4m_indep_gwas.txt
MYJOB_NAME=sld4m_indep
sbatch --export=NAME_LIST=$NAME_LIST,MYJOB_NAME=$MYJOB_NAME run_sld4m.sbatch
```

## T1 IDPs

```
# mkdir -p idp_lists
# cat ../../misc_data/supplementary_materials/supp_table_1.tsv |grep T1 | cut -f 1 | awk '{print "IDP-"$1}' > idp_lists/t1_list
# for i in `seq 1 5`; do echo "PC-"$i >> idp_lists/t1_list; done
# cd idp_lists
# split -l 50 -d t1_list t1_
# rm t1_list
# cd ../
for i in `ls idp_lists/t1_*`
do
  NAME_LIST=`pwd -P`/$i
  MYJOB_NAME=sld4m_t1
  echo sbatch --export=NAME_LIST=$NAME_LIST,MYJOB_NAME=$MYJOB_NAME run_sld4m.sbatch
done
```

## dMRI IDPs

```
# mkdir -p idp_lists
# cat ../../misc_data/supplementary_materials/supp_table_1.tsv |grep dMRI | cut -f 1 | awk '{print "IDP-"$1}' > idp_lists/dmri_list
# for i in `seq 1 9`; do echo "PC-"$i >> idp_lists/dmri_list; done
# cd idp_lists
# split -l 50 -d dmri_list dmri_
# rm dmri_list 
# cd ../
for i in `ls idp_lists/dmri_*`
do
  NAME_LIST=`pwd -P`/$i
  MYJOB_NAME=sld4m_dmri
  echo sbatch --export=NAME_LIST=$NAME_LIST,MYJOB_NAME=$MYJOB_NAME run_sld4m.sbatch
done
```

## Third round

```
# mkdir -p idp_lists_3rd
# ls /project2/haky/Data/BrainXcan/idp_gwas_3rd/trans_qtl.third_round_t1.chr22/* | sed 's#/project2/haky/Data/BrainXcan/idp_gwas_3rd/trans_qtl.third_round_t1.chr22/##g' | sed 's#.parquet##g' > idp_lists_3rd/t1_list
# cd idp_lists_3rd
# split -l 50 -d t1_list t1_
# rm t1_list
# cd ../
for i in `ls idp_lists_3rd/t1_*`
do
  NAME_LIST=`pwd -P`/$i
  MYJOB_NAME=sld4m_t1_3rd
  echo sbatch --export=NAME_LIST=$NAME_LIST,MYJOB_NAME=$MYJOB_NAME run_sld4m.sbatch
done
```

```
# mkdir -p idp_lists_3rd
# ls /project2/haky/Data/BrainXcan/idp_gwas_3rd/trans_qtl.third_round_dmri.chr22/* | sed 's#/project2/haky/Data/BrainXcan/idp_gwas_3rd/trans_qtl.third_round_dmri.chr22/##g' | sed 's#.parquet##g' > idp_lists_3rd/dmri_list
# cd idp_lists_3rd
# split -l 50 -d dmri_list dmri_
# rm dmri_list
# cd ../
for i in `ls idp_lists_3rd/dmri_*`
do
  NAME_LIST=`pwd -P`/$i
  MYJOB_NAME=sld4m_dmri_3rd
  echo sbatch --export=NAME_LIST=$NAME_LIST,MYJOB_NAME=$MYJOB_NAME run_sld4m.sbatch
done
```
