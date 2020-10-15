Here we perform ImageXcan analysis using observed phenotype and predicted IDPs.

* Input:
    - Observed phenotype table
    - Covariates
    - Predicted IDPs
    - Inclusion list

* Output
    - Summary statistic bhat, se, and z-score of the ImageXcan association test

* Association test:
    - Logistic regression (for binary trait). We **copy and modify** the solver implemented at [link](https://github.com/liangyy/haplotype-po/blob/master/scripts/logistic_gpu/logistic_gpu.py)
    - Linear regression (for quantitative trait). We implement the linear regression solver at `linear_regression_gpu.py`
