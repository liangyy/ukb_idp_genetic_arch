# prepare IDP table for Andrew Brown
# setwd('misc_data/download_some_matching_files/')
library(dplyr)
df = readRDS('cleanup_annot_our_idps.rds')
df = df[, -(1:6)]
df = df %>% select(FieldID, Field, Link, modality, ValueType, Units, Notes) %>% rename(Type = modality)
write.csv(df, 'idp_meta_table.csv')
