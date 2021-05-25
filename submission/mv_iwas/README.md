Here we run scripts in [MVIWAS](https://github.com/kathalexknuts/MVIWAS)

* 1000G reference panel 
* Download UKB IDP GWAS
* Prepare UKB IDP weights

Steps:

1. `Rscript prep_1000g_eur_list.R`
2. `bash submit_prep_1000g.sh` (note that for RCC, we need to separate download and other process ..)
3. `bash submit_ukb_idp_gwas.sbatch` (also, we need to separate download and others on RCC ..)
