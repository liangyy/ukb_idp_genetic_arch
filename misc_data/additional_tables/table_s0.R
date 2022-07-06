# table_s0
outdir <- '/Users/yanyuluchicago/Desktop/tmp/ukb_idp/additional_tables'
df <- read.csv('~/Desktop/tmp/ukb_idp/additional_tables/demographic.csv')
ff <- df$sex[df$X == 'sum']
mm <- df$sex[df$X == 'count_zeros']
ff0 <- ff / (ff + mm)
mm0 <- mm / (ff + mm)
d2 <- data.frame(
  variable = c('age', 'sex'),
  information = c(
    paste0(signif(df$age[df$X == 'mean'], digits = 3), ' (', signif(df$age[df$X == 'std'], digits = 3), ')'),
    paste0(signif(mm0 * 100, digits = 3), '%/', signif(ff0 * 100, digits = 3), '%')
  ),
  detail = c(
    paste0('range: ', df$age[df$X == 'min'], '--', df$age[df$X == 'max']), 
    paste0('raw count: ', mm, '/', ff)
  )
)
bold <- function(x) {paste('{\\textbf{',x,'}}', sep ='')}
print(
  xtable::xtable(d2, display=c("s", "s", "s", "s"), caption = '\\textbf{Demographic information of the IDP cohort.} A summary of age and sex of the selected 24,409 individuals (IDP cohort) is shown. For age, \\textbf{information} shows the average age with standard deviation of age in parentheses and \\textbf{detail} shows the range of the age (min and max). For sex, \\textbf{information} shows the fraction of each age among the cohort and \\textbf{detail} shows the raw counts (female\'s is listed first and followed by male\'s)', label = 'S-tab:ukb_cohort', digits = 10), 
  math.style.exponents = TRUE, 
  include.rownames = FALSE, 
  file = paste0(outdir, '/table_s0.tex'), 
  sanitize.colnames.function = bold,
  booktabs = TRUE
)



