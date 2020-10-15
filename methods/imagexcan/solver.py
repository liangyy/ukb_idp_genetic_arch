import numpy as np
import logistic_regression_gpu
import linear_regression_gpu

def logistic_regression(y, X, C, device=None):
    if not check_binary(y):
        raise ValueError('Input y is not binary. Cannot do logistic regression.')
    solver = logistic_regression_gpu.BatchLogisticSolver()
    bhat, se, success = solver.batchIRLS(
        X, y, C, 
        device=device, use_mask=True, min_prob=1e-20
    )
    pval, _, bhat, _ = solver.stat_test(bhat, se, success)
    return bhat, pval

def linear_regression(y, X, C, device=None):
    solver = linear_regression_gpu.LinearRegSolver()
    bhat, se, dof = solver.solve(
        X, y, C, 
        device=device
    )
    pval, _ = solver.stat_test(bhat, se, dof)
    return bhat, pval
