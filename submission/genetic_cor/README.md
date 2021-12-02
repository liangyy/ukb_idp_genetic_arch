Lists:

* GTEx-GWAS list: `/gpfs/data/im-lab/nas40t2/yanyul/GitHub/ukb_idp_genetic_arch/submission/simagexcan/gtex_gwas_list.txt`
* Psychiatric list (full): `/gpfs/data/im-lab/nas40t2/yanyul/GitHub/ukb_idp_genetic_arch/misc_data/preprocess_psychiatric_traits/trait_list.txt`
* Psychiatric list (no PD_Nalls_2019, needed for the formatting step): `psychiatric_traits_no_PD_Nalls_2019.txt`
* Psychiatric list (PD_Nalls_2019, needed for the formatting step): `psychiatric_traits_PD_Nalls_2019.txt`

3rd round submission.

Example

```
# skip formatting GWAS since we could re-use the formatted GWAS from 2nd round.
# ARGS1: pheno list
# ARGS2: gwas name tag
# ARGS3: IDP tag
# ARGS4: IDP list

# gwas name tag: psychiatric_3rd (use with all psychiatric traits), gtex_gwas_3rd 
# idp tags: t1_3rd, dmri_3rd
# no need to use no_PD_Nalls since it is only for formatting

bash submit_run.sh \
  /gpfs/data/im-lab/nas40t2/yanyul/GitHub/ukb_idp_genetic_arch/submission/simagexcan/gtex_gwas_list.txt \
  gtex_gwas_3rd \
  t1_3rd \
  /gpfs/data/im-lab/nas40t2/yanyul/GitHub/ukb_idp_genetic_arch/submission/genetic_cor/third_round_t1.txt
```

Fourth round (residual only)

```
bash submit_run.sh \
  /gpfs/data/im-lab/nas40t2/yanyul/GitHub/ukb_idp_genetic_arch/submission/simagexcan/gtex_gwas_list.txt \
  gtex_gwas_4th \
  t1_4th \
  /gpfs/data/im-lab/nas40t2/yanyul/GitHub/ukb_idp_genetic_arch/submission/genetic_cor/fourth_round_t1_residual.txt

bash submit_run.sh \
  /gpfs/data/im-lab/nas40t2/yanyul/GitHub/ukb_idp_genetic_arch/misc_data/preprocess_psychiatric_traits/trait_list.txt \
  psychiatric_4th \
  t1_4th \
  /gpfs/data/im-lab/nas40t2/yanyul/GitHub/ukb_idp_genetic_arch/submission/genetic_cor/fourth_round_t1_residual.txt
```

Test permutation

```
pheno=gtex_gwas_4th_w_perm
gwas=pgc.scz2
qsub -v NAME=$pheno,GWASNAME=$gwas format_gwas.qsub

bash submit_run_w_perm.sh \
  $gwas \
  $pheno \
  t1 \
  /gpfs/data/im-lab/nas40t2/yanyul/GitHub/ukb_idp_genetic_arch/submission/genetic_cor/fourth_round_t1_residual.txt \
  10 1 \
  /gpfs/data/im-lab/nas40t2/yanyul/data/ldblock/fourier_ls-all.bed
```
