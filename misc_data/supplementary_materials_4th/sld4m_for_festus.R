# SLD4M for Festus

backup_dir = '~/Documents/repo_backup/ukb_idp_genetic_arch/misc_data/supplementary_materials_4th/'

setwd(backup_dir)
library(ggplot2)
library(dplyr)
theme_set(theme_bw(base_size = 15))
source('https://gist.githubusercontent.com/liangyy/43912b3ecab5d10c89f9d4b2669871c9/raw/3ca651cfa53ffccb8422f432561138a46e93710f/my_ggplot_theme.R')

options(stringsAsFactors = F)

source('rlib.R')

foldern = '~/Documents/repo/GitHub/ukb_idp_genetic_arch/misc_data/supplementary_materials_4th/polygenicity'
dir.create(foldern)

ind = read.table('../ld4m/ld4m_indep_gwas.txt')$V1
dd3 = list()
url <- 'https://data.broadinstitute.org/alkesgroup/LDSCORE/independent_sumstats'
for(i in ind) {
  fn = paste0('~/Desktop/tmp/ukb_idp/ld4m/indep_gwas/', i, '.sld4m_all.csv')
  if(!file.exists(fn)) {
    next
  }
  tmp = read.csv(fn)
  dd3[[length(dd3) + 1]] = tmp %>% mutate(phenotype = i, type = 'independent', url = paste0(url, '/', i, '.sumstats.gz'))
}
ext = read.table('../ld4m/ld4m_external_traits.txt')$V1
url <- 'https://data.broadinstitute.org/alkesgroup/UKBB'
for(i in ext) {
  fn = paste0('~/Desktop/tmp/ukb_idp/ld4m/external_gwas/', i, '.sld4m_all.csv')
  if(!file.exists(fn)) {
    next
  }
  tmp = read.csv(fn)
  dd3[[length(dd3) + 1]] = tmp %>% mutate(phenotype = i, type = 'external', url = paste0(url, '/', i, '.sumstats.gz'))
}
dd3 = do.call(rbind, dd3)

ss = read.table('../ld4m/sample_size.txt', header = T)
dd3_sub = dd3 %>% filter(Var1 == 'Manual_aggregated')

df_poly_other = left_join(dd3_sub, ss, by = c('phenotype' = 'trait')) 
write.csv(df_poly_other, paste0(foldern, '/sld4m.yanyul_reproduce.csv'), row.names = FALSE)
