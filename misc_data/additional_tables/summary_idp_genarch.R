library(dplyr)

outdir <- '/Users/yanyuluchicago/Desktop/tmp/ukb_idp/additional_tables'
input <- '/Users/yanyuluchicago/Documents/repo_backup/ukb_idp_genetic_arch/misc_data/supplementary_materials_4th/result_tables/supp_genetic_arch.csv'
input_r2 <- '/Users/yanyuluchicago/Documents/repo_backup/ukb_idp_genetic_arch/misc_data/supplementary_materials_4th/supp_table_2.tsv'
df <- read.csv(input)
df_r2 <- read.delim2(input_r2) %>% filter(model_name == 'ridge')
df_r2$R2 <- as.numeric(df_r2$R2)
df <- left_join(df, df_r2, by = 'IDP')
head(df)
labels <- table(df$subtype) %>% as.data.frame() %>% filter(Freq >= 5) %>% pull(Var1)

# summarize each subtype by 
# T1/dMRI, subtype, 
# heritability of PC, mean heritability, sd of heritability, 
# Me of PC, mean Me, sd of Me, 
# R2(PC), mean R2, sd R2
res <- list()
for(label in labels) {
  sub <- df %>% filter(subtype == label)
  pc <- sub[ substr(sub$IDP, 1, 2) == 'PC', ]
  if(nrow(pc) == 0) 
    res0 <- data.frame(modality = sub$IDP_type[1], subtype = sub$subtype[1], h2 = NA, Me = NA, R2 = NA)
  else 
    res0 <- pc %>% select(IDP_type, subtype, h2, Me, R2)
  colnames(res0) <- c('modality', 'subtype', 'h2(PC)', 'Me(PC)', 'R2(PC)')
  res1 <- sub[ substr(sub$IDP, 1, 2) != 'PC', ] %>% 
    summarize(
      `h2(mean)` = mean(h2), 
      `Me(mean)` = mean(Me),
      `R2(mean)` = mean(R2),
      `h2(s.d.)` = sd(h2), 
      `Me(s.d.)` = sd(Me),
      `R2(s.d.)` = sd(R2),
      `Number of IDPs` = n())
  res[[length(res) + 1]] <- cbind(res0, res1)
}
res <- do.call(rbind, res)
n <- 3
res <- res %>% 
  mutate(
    h2 = paste0(signif(`h2(mean)`, digits = n), ' (', signif(`h2(s.d.)`, digits = n), ')'),
    Me = paste0(formatC(`Me(mean)`, digits = n - 1, format = 'e'), ' (', formatC(`Me(s.d.)`, digits = n - 1, format = 'e'), ')'),
    R2 = paste0(signif(`R2(mean)`, digits = n), ' (', signif(`R2(s.d.)`, digits = n), ')'),
    `h2(PC)` = signif(`h2(PC)`, digits = n),
    `Me(PC)` = formatC(`Me(PC)`, digits = n - 1, format = 'e'),
    `R2(PC)` = signif(`R2(PC)`, digits = n)) %>%
  select(
    modality, subtype, `Number of IDPs`,
    `h2(PC)`, h2,
    `Me(PC)`, Me,
    `R2(PC)`, R2) %>% 
  arrange(desc(modality), subtype) # %>% 
  # filter(substr(subtype, 1, 1) != 'w')
write.csv(res, paste0(outdir, '/summary_idp_genarch.csv'), row.names = FALSE)
