# table 1
library(dplyr)

outdir <- '/Users/yanyuluchicago/Desktop/tmp/ukb_idp/additional_tables'
dd <- read.csv(paste0(outdir, '/summary_idp_genarch.csv'))
colnames(dd)[3] <- 'Number of IDPs'
colnames(dd) <- stringr::str_replace(colnames(dd), '.PC.', '(PC)')
print(
  xtable::xtable(dd, display=rep("s", 10), caption = '\\textbf{Demographic information of the IDP cohort.} A summary of age and sex of the selected 24,409 individuals (IDP cohort) is shown. For age, \\textbf{information} shows the average age with standard deviation of age in parentheses and \\textbf{detail} shows the range of the age (min and max). For sex, \\textbf{information} shows the fraction of each age among the cohort and \\textbf{detail} shows the raw counts (female\'s is listed first and followed by male\'s)', label = 'S-tab:ukb_cohort', digits = 10), 
  math.style.exponents = TRUE, 
  include.rownames = FALSE, 
  file = paste0(outdir, '/table_1.tex'), 
  sanitize.colnames.function = bold,
  booktabs = TRUE
)
