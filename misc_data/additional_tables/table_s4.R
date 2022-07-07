# table_s4
library(dplyr)

outdir <- '/Users/yanyuluchicago/Desktop/tmp/ukb_idp/additional_tables'
df <- read.delim2('/Users/yanyuluchicago/Documents/repo_backup/ukb_idp_genetic_arch/misc_data/supplementary_materials_4th/supp_table_3.tsv')
write.csv(df, paste0(outdir, '/table_s4.csv'), row.names = FALSE)