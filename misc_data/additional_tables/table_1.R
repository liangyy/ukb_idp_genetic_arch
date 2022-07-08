# table 1
library(dplyr)

outdir <- '/Users/yanyuluchicago/Desktop/tmp/ukb_idp/additional_tables'
dd <- read.csv(paste0(outdir, '/summary_idp_genarch.csv'))
colnames(dd)[3] <- 'Number of IDPs'
colnames(dd) <- stringr::str_replace(colnames(dd), '.PC.', '(PC)')

bold <- function(x) {paste('{\\textbf{',x,'}}', sep ='')}
print(
  xtable::xtable(dd, display=rep("s", 10), caption = '\\textbf{A summary of the genetic architecture of prediction performance of brain IDPs.} 
  This table summarizes the genetic architecture estimates, heritability and polygenicity (in terms of $M_e$ \\citep{oconnor:2019}), and the prediction performance of ridge models (in terms of $R^2$) of brain IDPs organized by IDP subtypes (defined in table~\\ref{S-tab:idp_table}). 
  \textbf{modality} and \\textbf{subtype} show the modality of subtype name (see definitions in table~\ref{S-tab:idp_table}). 
  \\textbf{Number of IDPs} shows the number of brain IDPs falling in the subtype. 
  \\textbf{h2}, \\textbf{Me}, and \\textbf{R2} show the average heritability, polygenicity, and prediction performance among the IDPs (excluding principal component) within the subtype and the corresponding standard deviation is shown inside the parentheses.
  \\textbf{h2(PC)}, \\textbf{Me(PC)}, and \\textbf{R2(PC)} show the heritability, polygenicity, and prediction performance of the principal component of the subtype. ', label = 'S-tab:summary_of_gen_arch', digits = 10), 
  math.style.exponents = TRUE, 
  include.rownames = FALSE, 
  file = paste0(outdir, '/table_1.tex'), 
  sanitize.colnames.function = bold,
  booktabs = TRUE
)
