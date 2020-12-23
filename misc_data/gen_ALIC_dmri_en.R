library(dplyr)

df_d = arrow::read_parquet('~/Desktop/tmp/ukb_idp/gw_elastic_net_dmri/gw_lasso_beta.parquet')
# df_t = arrow::read_parquet('~/Desktop/tmp/ukb_idp/gw_elastic_net_t1/gw_lasso_beta.parquet')
annot = readRDS('misc_data/download_some_matching_files/annot_dmri_idps.rds')
idps = annot[ annot$position == 'anterior limb of internal capsule', ]
cols = c('snpid', 'a0', 'a1', 'chr', paste0('IDP-', idps$FieldID))
df_d = df_d[, cols]
df_d = df_d[ rowSums(df_d[ -(1:4)]) != 0,  ]
df_annot = inner_join(
  data.frame(idp = colnames(df_d)[-(1:4)]), 
  annot %>% mutate(idp = paste0('IDP-', FieldID)) %>% select(idp, Field), 
  by = 'idp'
)
gz1 = gzfile('misc_data/weights.ALIC.dmri_en.tsv.gz', "w")
write.table(df_d, gz1, quo = F, col = T, row = F, sep = '\t')
close(gz1)
gz1 = gzfile('misc_data/annot.ALIC.dmri_en.tsv.gz', "w")
write.table(df_annot, gz1, quo = F, col = T, row = F, sep = '\t')
close(gz1)

df_snp = list()
for(cc in cols) {
  tmp = df_d[, c('snpid', 'a0', 'a1', 'chr', cc)] 
  tmp = tmp[ tmp[, cc] != 0, ] 
  df_snp[[length(df_snp) + 1]] = tmp[, c('snpid', 'a0', 'a1', 'chr')]
}
df_snp = do.call(rbind, df_snp)

gz1 = gzfile('misc_data/snp_list.ALIC.dmri_en.tsv.gz', "w")
write.table(df_snp, gz1, quo = F, col = T, row = F, sep = '\t')
close(gz1)

mr = readRDS('~/Desktop/tmp/ukb_idp/mr_summary_121520/MR.gtex_gwas_dmri.IDP-25361_x_pgc.scz2.rds')
mr$idp2pheno$data
df_snp2 = rbind(
  data.frame(snpid = mr$idp2pheno$data$SNP, direction = 'idp2pheno'),
  data.frame(snpid = mr$pheno2idp$data$SNP, direction = 'pheno2idp')
)
gz1 = gzfile('misc_data/mr_snp_list.ALIC.dmri_en.tsv.gz', "w")
write.table(df_snp2, gz1, quo = F, col = T, row = F, sep = '\t')
close(gz1)