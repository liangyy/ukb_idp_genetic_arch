# setwd('misc_data/supplementary_materials')

# here we form the list of GWAS for S-BrainXcan analysis

library(dplyr)

# gtex-gwas
df = read.csv('~/Downloads/gtex-gwas-full - Sheet1.csv')
other_pheno = c('BMI_EUR', 'BMI_Active_Inds', 'Height', 'Standing_Height_UKB', 'BMI_UKB')
exclude_psych = c('ASD_PGC_2017', 'ADHD_PGC_2017', 'BMI_Locke_2018')
df_selected = df %>% filter(Category == 'Psychiatric-neurologic' | new_Phenotype %in% other_pheno) %>% filter(Cases > 1000 | is.na(Cases))


# psychiatric
df = read.csv('~/Downloads/psychiatric_traits_downloaded - Sheet1.csv')

out = rbind(
  df_selected %>% select(Phenotype, Tag, new_abbreviation, Sample_Size, Portal, Pheno_File) %>% 
    rename(phenotype = Phenotype, phenotype_id = Tag, short_name = new_abbreviation, sample_size = Sample_Size, portal = Portal, filename = Pheno_File),
  df %>% 
    mutate(sample_size = c(53293, 412144, 46351, 24409, 23261, 686648, 766345, 44960, 258327, 106764, 383873, 455284, 69826, 71881), short_name = trait_id) %>%
    select(trait, trait_id, short_name, sample_size, filename, source) %>% 
    rename(phenotype = trait, phenotype_id = trait_id, portal = source)
)
out = out %>% filter(!phenotype_id %in% exclude_psych)
out$portal[out$portal == 'PGC'] = 'https://www.med.unc.edu/pgc/results-and-downloads'
out$portal[out$portal == 'http://biobank.ctsu.ox.ac.uk/'] = 'http://www.nealelab.is/uk-biobank'
write.table(out, 'supp_table_4.tsv', quote = F, row = F, col = T, sep = '\t')
