library(optparse)

option_list <- list(
    make_option(c("-g1", "--group1"), type="character", default=NULL,
                help="Output group1",
                metavar="character"),
    make_option(c("-g2", "--group2"), type="character", default=NULL,
                help="Output group2",
                metavar="character"),
    make_option(c("-i", "--input"), type="character", default=NULL,
                help="Input FAM file",
                metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

set.seed(2021)

samples <- data.table::fread(opt$input, data.table = FALSE, header = FALSE)$V1
n1 <- 10000
n2 <- length(samples) - n1
g1_idx <- sample(1 : length(samples), n1)
g1_ind <- (1 : length(samples)) %in% g1_idx
g1 <- samples[g1_ind]
g2 <- samples[g2_ind]
message(glue::glue('Number of samples in group1 = {length(g1)}'))
message(glue::glue('Number of samples in group2 = {length(g2)}'))

write.table(
  data.frame(x = g1, y = g1),
  opt$group1,
  col = FALSE, row = FALSE, sep = '\t')
write.table(
  data.frame(x = g2, y = g2),
  opt$group2,
  col = FALSE, row = FALSE, sep = '\t')
