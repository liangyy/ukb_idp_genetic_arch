options(stringsAsFactors = F)
dd = c(read.table('ld4m_external_traits.txt')$V1, read.table('ld4m_indep_gwas.txt')$V1)
df_n = list()
for(i in dd) {
  message(i)
  tmp = data.table::fread(cmd = paste0('zcat /project2/haky/yanyul/data/SLD4M/external_gwas/', i, '.sumstats.gz'))
  df_n[[length(df_n) + 1]] = data.frame(trait = i, mean_sample_size = mean(tmp$N, na.rm = T))
}
df_n = do.call(rbind, df_n)

write.table(df_n, 'sample_size.txt', quo = F, row = F, col = T, sep = '\t')
