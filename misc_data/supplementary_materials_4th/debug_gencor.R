library(dplyr)
options(stringsAsFactors=F)
source('rlib.R')
cmd = "cat ~/scratch/ukb_idp/mv_iwas/ukb_idp_gwas/25006.tmp2.txt|awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9}'"; idp18 = fread(cmd = cmd, data.table = F)
meta = fread('~/labshare/ukb_idp/subset_genotypes/IDP_HM3_finalPheno.chr22.bim', data.table = F)
my = arrow::read_parquet('~/scratch/ukb_idp/idp_gwas_4th/trans_qtl.fourth_round.t1_w_pc.chr22/IDP-25006.parquet')
idp1822 = idp18 %>% filter(CHR == 22)
m1 = inner_join(idp1822, my, by = c('SNP' = 'variant_id'), suffix = c('.18', '.my'))

# The alleles are defined such that a1 is the reference allele, and a2 is the alternative allele; hence, a2 will often be the minor allele, but not always. The effect allele (in the GWAS linear model) is a2, meaning that the sign of regression beta relates to the a2 allele count. AF relates to the frequency of a2. Hence MAF (minor allele frequency) = min(AF,1-AF). The reference human assembly is GRCH37/hg19. Our rsids are mostly covered by dbSNPs 147.

kk = fread(cmd = 'zcat ~/labshare/tmp/ukb_idp/idp_gwas_2021/0005.txt.gz', data.table = F)
kk22 = kk %>% filter(chr == 22)
m2 = inner_join(my, kk22, by = c('variant_id' = 'rsid'))

tmp = read.table('~/Desktop/tmp/ukb_idp/idp_gwas_2021/temp.chr22.idp25006-25005.txt')
tmp2 = read.table('~/Desktop/tmp/ukb_idp/idp_gwas_2021/temp.chr22.idp25006-25006.txt')
head(tmp)
plot(-log10(tmp$pval), tmp$pval..log10.)
plot(-log10(tmp2$pval), tmp2$pval..log10.)


library(ggplot2)
tmp = tmp %>% mutate(direction = (V5 == a2))
tmp2 = tmp2 %>% mutate(direction = (V5 == a2))
tmp %>% ggplot() + geom_point(aes(x = b, y = beta)) + facet_wrap(~direction)
tmp2 %>% ggplot() + geom_point(aes(x = b, y = beta)) + facet_wrap(~direction)

# library(dplyr)
# options(stringsAsFactors=F)
# tmp = list()
# for(i in 1 : 22) {
#   t = arrow::read_parquet(paste0('~/scratch/ukb_idp/idp_gwas_4th/trans_qtl.fourth_round.t1_w_pc.chr', i, '/IDP-25006.parquet'))
#   tmp[[length(tmp) + 1]] = t
# }
# tmp = do.call(rbind, tmp)
# meta = fread('~/labshare/ukb_idp/subset_genotypes/IDP_HM3_finalPheno.merged_all.bim', data.table = F)
# colnames(meta) = c('chr', 'variant_id', 'nn', 'pos', 'alt', 'ref')
# tmp = left_join(tmp, meta %>% select(variant_id, alt, ref), by = 'variant_id')
# tmp %>% write.table('my.IDP-25006', quo = F, row = F, col = T, sep = '\t')
library(data.table)
d1 = fread('zcat < ~/Desktop/tmp/ukb_idp/idp_gwas_2021/formatted.0006.sumstats.gz', data.table = F)
d2 = fread('zcat < ~/Desktop/tmp/ukb_idp/idp_gwas_2021/formatted.25006.sumstats.gz', data.table = F)
mm = inner_join(d1, d2, by = 'SNP', suffix = c('.0006', '.25006'))
mm = mm %>% mutate(direction = A1.0006 == A1.25006)
mm[sample(1 : nrow(mm), 10000, replace = F), ] %>% ggplot() + geom_point(aes(x = Z.0006, y = Z.25006)) + facet_wrap(~direction)

ii = load_idp_annot()
kk2 = readRDS('s_bxcan/dataframe_full.gencor.rds') %>% filter(phenotype == 'UKB_20016_Fluid_intelligence_score') %>% inner_join(ii, by = 'IDP')
kk = readRDS('s_bxcan/dataframe_full.gencor.rds') %>% filter(phenotype == 'Intelligence_CTG') %>% inner_join(ii, by = 'IDP')
kk3 = readRDS('s_bxcan/dataframe_full.sbxcan.rds') %>% filter(phenotype == 'UKB_20016_Fluid_intelligence_score' | phenotype == 'Intelligence_CTG') %>% inner_join(ii, by = 'IDP')
