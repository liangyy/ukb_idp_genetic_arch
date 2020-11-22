# grab the related GWASs of these phenotypes of interest
# from open GWAS database (https://gwas.mrcieu.ac.uk/)

library(dplyr)

pheno_interest = c('weekly_alcohol', 'recurrent_depressive_disorder', 'parent_depression', 'parent_AD', 'handedness', 'daily_coffee', 'daily_cigarettes')

res = list()

# weekly_alcohol
res[['weekly_alcohol']] = c('ukb-b-5779', 'ukb-a-25')

# recurrent_depressive_disorder
res[['recurrent_depressive_disorder']] = c('ukb-b-12064', 'ebi-a-GCST005902', 'ukb-d-20544_11')

# parent_depression
res[['parent_depression']] = c('ukb-b-12064', 'ebi-a-GCST005902', 'ukb-d-20544_11')

# parent_AD
res[['parent_AD']] = c('ebi-a-GCST005921', 'ieu-b-2', 'ukb-b-14699', 'ukb-b-323')

# handedness
res[['handedness']] = c('ukb-d-1707_1', 'ukb-a-36')


# daily_coffee
res[['daily_coffee']] = c('ukb-b-5237', 'ukb-b-9508')

# daily_cigarettes
res[['daily_cigarettes']] = c('ieu-b-25', 'ukb-b-6019', 'ieu-a-961', 'ukb-b-469', 'ukb-a-342', 'ukb-b-20261', 'ukb-a-236')

# form data frame
oo = list()
for(n in names(res)) {
  oo[[length(oo) + 1]] = data.frame(phenotype = n, gwas_code = res[[n]])
}
oo = do.call(rbind, oo)

# add meta information
ao <- TwoSampleMR::available_outcomes()
oo = left_join(oo, ao, by = c('gwas_code' = 'id'))

# assign gwas id: phenotype + 1 : n
oo = oo %>% group_by(phenotype) %>% mutate(id = paste0(phenotype, '_', 1 : n()))

# save
write.csv(oo, 'pheno_of_interest.meta_table.csv')
