#! /usr/bin/env Rscript

library(arrow)
library(data.table)
library(tidyverse)
out.dir <- "/Volumes/im-lab/nas40t2/festus/brainxcan"


# dmri
idps <- list.files(path = glue::glue("{out.dir}/results.all/dmri"),
                    pattern = ".*.txt", full.names = T)
weights = ""
for (idp in idps) {
  
  gene <- str_remove(basename(tools::file_path_sans_ext(idp))," ")
  
  betas <- fread(idp) 
  if(nrow(betas) == 0) next
  betas <- betas %>% 
    setnames(names(.), c("chr","snpid", "pos", "a1", "a0", gene)) %>% 
    dplyr::select(snpid,a0,a1,chr,all_of(gene))
  
  ## bind them together
  if (exists('weights') && is.data.frame(get('weights'))) {
    weights <- weights %>% 
      inner_join(betas, by = c("snpid","a0","a1","chr"))
  } else {
    weights <- betas
  }
}
write_parquet(weights, glue::glue("{out.dir}/models/fourth_round.dmri_w_pc.parquet"))

# t1
idps <- list.files(path = glue::glue("{out.dir}/results.all/t1"),
                    pattern = ".*.txt", full.names = T)
weights = ""
for (idp in idps) {
  
  gene <- str_remove(basename(tools::file_path_sans_ext(idp))," ")
  
  betas <- fread(idp) 
  if(nrow(betas) == 0) next
  betas <- betas %>% 
    setnames(names(.), c("chr","snpid", "pos", "a1", "a0", gene)) %>% 
    dplyr::select(snpid,a0,a1,chr,all_of(gene))
  
  ## bind them together
  if (exists('weights') && is.data.frame(get('weights'))) {
    weights <- weights %>% 
      inner_join(betas, by = c("snpid","a0","a1","chr"))
  } else {
    weights <- betas
  }
}
write_parquet(weights, glue::glue("{out.dir}/models/fourth_round.t1_w_pc.parquet"))
