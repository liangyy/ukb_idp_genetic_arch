# generate map between selected phenotypes and open GWAS data.
# see the definition of the selected phenotypes below
pheno_interest = c('weekly_alcohol', 'recurrent_depressive_disorder', 'parent_depression', 'parent_AD', 'handedness', 'daily_coffee', 'daily_cigarettes')
pheno_bcc = c('wbc', 'rbc', 'platelet', 'lymphocyte', 'monocyte', 'neutrophil', 'eosinophil', 'basophil')
pheno_imgx = c('height', 'bmi', pheno_interest, pheno_bcc)

gwas_code_interest = c('ukb-b-5779', 'ukb-b-12064', 'ukb-b-12064', 'ebi-a-GCST005921', 'ukb-a-36', 'ukb-b-5237', 'ieu-b-25')
gwas_code_bcc = c('ieu-b-30', 'ebi-a-GCST004601', 'ukb-d-30080_irnt', 'ebi-a-GCST004627', 'ukb-d-30130_irnt', 'ebi-a-GCST004629', 'ebi-a-GCST004606', 'ieu-b-29')
gwas_code = c('ukb-b-10787', 'ukb-b-19953', gwas_code_interest, gwas_code_bcc)
df_ogwas2 = data.frame(phenotype = pheno_imgx, gwas_code = gwas_code)
saveRDS(df_ogwas2, 'misc_data/selected_pheno_to_open_gwas.rds')