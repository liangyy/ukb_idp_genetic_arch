standardize = function(x) {
  apply(x, 2, function(y) { (y - mean(y)) / sd(y) })
}

corr_cat_cont = function(categ, cont, ...) {
  mod = lmer(cont ~ (1 | categ), ...)
  if(!isSingular(mod, tol = 1e-05)) {
    test = lmerTest::ranova(mod)
    return(test$`Pr(>Chisq)`[2])
  } else {
    return(1)
  }
}
