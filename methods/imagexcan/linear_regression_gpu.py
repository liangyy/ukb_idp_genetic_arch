import scipy.stats
import torch
import statsmodels.api as sm
import numpy as np

class LinearRegSolver:
    def __init__(self):
        self.info = 'Fit y ~ x + covariates for every x_i ~ y pair.'
    
    @staticmethod
    def solve(X, y, C, device=None):
        q, _ = torch.qr(C)
        y_ = y - torch.matmul(q, torch.matmul(q.T, y))
        X_ = X - torch.matmul(q, torch.matmul(q.T, X))
        xty = torch.matmul(X_.T, y_)
        xtx = torch.pow(X_, 2).sum(axis=0)
        bhat = xty / xtx
        dof = X.shape[0] - C.shape[1] - 1
        sigma2 = torch.pow(torch.unsqueeze(y_, 1) - torch.einsum('ij,j->ij', X_, bhat), 2).sum(axis=0) / dof 
        se = torch.sqrt(sigma2 / xtx)
        
        bhat = bhat.numpy()
        se = se.numpy()
        
        return bhat, se, dof
    
    @staticmethod    
    def stat_test(bhat, se, dof):
        tval = bhat / se
        pval = scipy.stats.t.sf(np.abs((bhat / se)), df=dof) * 2
        return pval, tval
    
def test_solver(y, x, covar, return_val=False):
    '''
    Input y, x, covar are numpy array
    y: n x 1
    x: n x k
    covar: n x p
    Return: k bhat's and se's by LinearRegSolver and by statsmodels.api.OLS
    '''
    print('---- Testing LinearRegSolver ----')
    # LinearRegSolver
    solver = LinearRegSolver()
    bhat, se, dof = solver.solve(torch.Tensor(x), torch.Tensor(y), torch.Tensor(covar))
    pval, stat = solver.stat_test(bhat, se, dof)
    
    
    # statsmodels
    K = x.shape[1]
    bhat2 = np.zeros((K))
    stat2 = np.zeros((K))
    pval2 = np.zeros((K))
    for k in range(K):
        mod = sm.OLS(y, np.concatenate((x[:, k:(k+1)], covar), axis=1))
        res = mod.fit()
        bhat2[k] = res.params[0]
        stat2[k] = res.tvalues[0]
        pval2[k] = res.pvalues[0]
    se2 = bhat2 / stat2
    print('bhat mean abs difference = {}'.format(np.abs((bhat - bhat2)).mean()))
    print('se mean abs difference = {}'.format(np.abs((se - se2)).mean()))
    print('tval mean abs difference = {}'.format(np.abs((stat - stat2)).mean()))
    print('log p mean abs difference = {}'.format(np.abs((np.log(pval) - np.log(pval2))).mean()))
    
    if return_val is True:
        return bhat, pval, bhat2, pval2
        
