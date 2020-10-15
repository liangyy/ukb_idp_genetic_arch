import numpy as np
from linear_regression_gpu import test_solver as test_solver_lm
from logistic_regression_gpu import test_solver as test_solver_logit

x = np.random.rand(100, 7)
covar = np.random.rand(100, 10)
y = x[:, 5] * 0.1 + covar.sum(axis=1) + np.random.rand(100) * 2

test_solver_lm(y, x, covar)
test_solver_lm(y, x, np.concatenate((np.ones((100, 1)), covar), axis=1))

ybin = np.zeros(100)
ybin[ y > y.mean() ] = 1

test_solver_logit(ybin, x, covar)
test_solver_logit(ybin, x, np.concatenate((np.ones((100, 1)), covar), axis=1))
