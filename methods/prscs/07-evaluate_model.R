#! /usr/bin/env Rscript

library(data.table)
library(tidyverse)
source("/Users/festo/Documents/HakyIm/mesa-metabolomics-analysis/code/compare_expression.R")
out.dir <- "/Volumes/im-lab/nas40t2/festus/brainxcan"

# load predicted and observed
#t1
observed <- read_parquet(glue::glue("{out.dir}/genotype/test/fourth_round.t1_w_pc.parquet")) %>% dplyr::rename(IID = individual) %>% mutate(FID = IID) %>% 
  relocate(FID, .before = IID)
predicted <- read_parquet(glue::glue("{out.dir}/models/predicted_t1_idps.parquet")) %>% dplyr::rename(IID = indiv) %>% mutate(FID = IID) %>% 
  relocate(FID, .before = IID)


t1_list <- fn_compare_expression_matrices(observed, predicted, nplots = 0)
t1_corr <- data.frame(IDP=t1_list[[4]], Spearman_prscs = t1_list[[3]])

#dmri
observed <- read_parquet(glue::glue("{out.dir}/genotype/test/fourth_round.dmri_w_pc.parquet")) %>% dplyr::rename(IID = individual) %>% mutate(FID = IID) %>% 
  relocate(FID, .before = IID)
predicted <- read_parquet(glue::glue("{out.dir}/models/predicted_dmri_idps.parquet")) %>% dplyr::rename(IID = indiv) %>% mutate(FID = IID) %>% 
  relocate(FID, .before = IID)


dmri_list <- fn_compare_expression_matrices(observed, predicted, nplots = 0)
dmri_corr <- data.frame(IDP=dmri_list[[4]], Spearman_prscs = dmri_list[[3]])

## save prediction corr
fwrite(dmri_corr, file = glue::glue("{out.dir}/prscs_prediction_perf-dmri.txt"),
       sep = "\t")
fwrite(t1_corr, file = glue::glue("{out.dir}/prscs_prediction_perf-t1.txt"),
       sep = "\t")
