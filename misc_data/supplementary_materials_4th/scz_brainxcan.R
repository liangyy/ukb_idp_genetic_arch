# setwd('misc_data/supplementary_materials_4th/')

library(dplyr)
source('rlib.R')


# quick number for SCZ PGC 2020 brainxcan run
options(stringsAsFactors = F)
mm = read.csv('~/Desktop/tmp/ukb_idp/brainxcan_pipeline/psychiatric/SCZ_PGC_2020.sbrainxcan.csv')
mm = remove_probtrack_idp(mm)
mm_top = mm[ 1 : 10, ]
mrs = list()
for(i in 1 : nrow(mm_top)) {
  s1 = mm_top$IDP[i]
  s2 = tolower(mm_top$modality[i])
  fn = paste0('/Users/yanyul/Desktop/tmp/ukb_idp/brainxcan_pipeline/psychiatric/SCZ_PGC_2020_files/MR_results/', s2, '_', s1, '.MR_sumstats.tsv')
  mrs[[length(mrs) + 1]] = read.delim2(fn) %>% filter(method == 'ACAT meta-analysis') %>% mutate(IDP = s1)
}
mrs = do.call(rbind, mrs)
mm_top = mm_top %>% left_join(mrs, by = 'IDP')
mm_top %>% group_by(direction) %>% summarize(ntotal = n(), nsig = sum(pval.y < 0.05))
mm_top %>% filter(pval.y < 0.05)
