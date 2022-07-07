# table s2
library(dplyr)
outdir <- '/Users/yanyuluchicago/Desktop/tmp/ukb_idp/additional_tables'
input <- '/Users/yanyuluchicago/Documents/repo_backup/ukb_idp_genetic_arch/misc_data/supplementary_materials_4th/result_tables/supp_genetic_arch.csv'
df <- read.csv(input)
write.csv(df %>% select(-subtype, -region, -side), paste0(outdir, '/table_s2.csv'), row.names = FALSE)
