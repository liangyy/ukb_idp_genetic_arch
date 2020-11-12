import numpy as np
import torch

from tqdm import tqdm
import logistic_regression_gpu
import linear_regression_gpu
from util import check_binary

import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import numpy2ri
numpy2ri.activate()
r_base = importr('base')
r_susieR = importr('susieR')

def logistic_regression(y, X, C, batch_size=5, device=None):
    if not check_binary(y):
        raise ValueError('Input y is not binary. Cannot do logistic regression.')
    solver = logistic_regression_gpu.BatchLogisticSolver()
    blist = []
    plist = []
    batch_size = min(batch_size, X.shape[1])
    nbatch = X.shape[1] // batch_size
    if nbatch * batch_size < X.shape[1]:
        nbatch += 1
    
    for i in tqdm(range(nbatch)):
        start = batch_size * i
        end = min(batch_size * (i + 1), X.shape[0])
        
        bhat, se, success, diff = solver.batchIRLS(
                torch.Tensor(X[:, start:end]), torch.Tensor(y[:, np.newaxis]), torch.Tensor(C), 
            device=device, use_mask=True, min_prob=1e-20, maxiter=100
        )
        
        pval, _, bhat, _ = solver.stat_test(bhat, se, success)
        blist.append(bhat)
        plist.append(pval)
    
    bhat = np.concatenate(blist, axis=0)
    pval = np.concatenate(plist, axis=0)
    return bhat, pval

def linear_regression(y, X, C, device=None):
    solver = linear_regression_gpu.LinearRegSolver()
    bhat, se, dof = solver.solve(
        torch.Tensor(X), torch.Tensor(y), torch.Tensor(C), 
        device=device
    )
    pval, _ = solver.stat_test(bhat, se, dof)
    return bhat, pval

def run_susie_wrapper(z, R):
    mod = r_susieR.susie_rss(z, R)
    df_ = r_base.summary(mod).rx2('vars')
    with localconverter(ro.default_converter + pandas2ri.converter):
        df_ = ro.conversion.rpy2py(df_)
    return list(df_.variable_prob), list(df_.cs)

