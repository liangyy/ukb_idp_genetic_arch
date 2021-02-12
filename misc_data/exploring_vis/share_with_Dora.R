library(dplyr)
kk = read.csv('~/Desktop/tmp/ukb_idp/simagexcan/results_psychiatric_2nd/dmri.original.all_covar.w_pc.gw_elastic_net_beta_x_Neuroticism_CTG_x_simagexcan.csv')
kk2 = read.csv('~/Desktop/tmp/ukb_idp/simagexcan/results_psychiatric_2nd/t1.scaled.all_covar.w_pc.gw_elastic_net_beta_x_Neuroticism_CTG_x_simagexcan.csv')
kk = rbind(kk, kk2)
idps = read.delim2('misc_data/supplementary_materials/supp_table_1.tsv')
kk = left_join(
  kk, 
  idps %>%
    mutate(IDP = paste0('IDP-', ukb_field)) %>%
    select(IDP, anatomy, left_or_right, measurement_type, dmri_measure, t1_anatomy_group, notes),
  by = 'IDP'
)
write.csv(kk, 'misc_data/exploring_vis/example_brainxcan_results_EN_Neuroticism.csv')
