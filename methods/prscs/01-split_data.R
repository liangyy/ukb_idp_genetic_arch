#! /usr/bin/env Rscript

library(arrow)
library(data.table)
library(tidyverse)
set.seed(918)
source("/Users/festo/Documents/HakyIm/mesa-metabolomics-analysis/code/compare_expression.R")
out.dir <- "/Volumes/im-lab/nas40t2/festus/brainxcan"

## Split the genotype
in.fam <- fread("/Volumes/im-lab/nas40t2/yanyul/ukb_idp/subset_genotypes/IDP_HM3_finalPheno.chr1.fam") %>% 
  setnames(., names(.), c("FID","IID","Father","Mother","Sex","Phenotype")) %>% 
  data.frame()

split1<- sample(c(rep(0, 0.8 * nrow(in.fam)), rep(1, 0.2 * nrow(in.fam))))
table(split1)
training_set <- in.fam[split1 == 0,] %>% select(FID,IID)
testing_set <- in.fam[split1 == 1,] %>% select(FID,IID)

fwrite(training_set, file = glue::glue("{out.dir}/train.ind.txt"), sep = " ",
       col.names = F)
fwrite(testing_set, file = glue::glue("{out.dir}/test.ind.txt"), sep = " ",
       col.names = F)

## split phenotype
#dmri
pq <- "/Volumes/im-lab/nas40t2/yanyul/ukb_idp/fourth_round_idp_preprocessing/fourth_round.dmri_w_pc.parquet"
dta = read_parquet(pq)

train.dta <- dta %>% 
  filter(individual %in% training_set$IID)
test.dta <- dta %>% 
  filter(individual %in% testing_set$IID)

write_parquet(train.dta, glue::glue("{out.dir}/genotype/train/fourth_round.dmri_w_pc.parquet"))

write_parquet(test.dta, glue::glue("{out.dir}/genotype/test/fourth_round.dmri_w_pc.parquet"))

#t1
pq <- "/Volumes/im-lab/nas40t2/yanyul/ukb_idp/fourth_round_idp_preprocessing/fourth_round.t1_w_pc.parquet"
dta = read_parquet(pq)

train.dta <- dta %>% 
  filter(individual %in% training_set$IID)
test.dta <- dta %>% 
  filter(individual %in% testing_set$IID)

write_parquet(train.dta, glue::glue("{out.dir}/genotype/train/fourth_round.t1_w_pc.parquet"))

write_parquet(test.dta, glue::glue("{out.dir}/genotype/test/fourth_round.t1_w_pc.parquet"))


