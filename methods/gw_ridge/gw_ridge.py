import numpy as np

# for debugging
import pdb

class blupRidgeSolver:
    
    def __init__(self, grm=None, y=None, theta_g_grid=None, inner_cv_fold=5):
        if theta_g_grid is None:
            self.theta_g_grid = np.array([0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95])
        self.inner_cv_fold = inner_cv_fold
        self.grm = grm
        self.y = y
        self.ny = y.shape[1]
        
    def get_partition(self, train_idx):
        idx = np.arange(self.inner_cv_fold)
        partition_pool = []
        size_part = train_idx.shape[0] // self.inner_cv_fold
        rest = train_idx.shape[0] - self.inner_cv_fold * size_part
        for i in range(self.inner_cv_fold):
            partition_pool += [i] * size_part
            if i < rest:
                partition_pool += [i]
        partition_pool = np.array(partition_pool)
        partitions = [] # each element will be a pair of training idx and testing idx
        for i in range(self.inner_cv_fold):
            part_test_idx = train_idx[partition_pool == i]
            part_train_idx = train_idx[partition_pool != i]
            partitions.append((part_train_idx, part_test_idx))
        return partitions
    
    def train(self, train_idx, test_idx, theta_g, subset_y_idx=None):
        '''
        return mse and predicted y
        formula:
            M = (1 - theta_g) * np.eye(train_idx.shape[0]) + theta_g * grm[train_idx, train_idx]
            ypred = theta_g * (grm[test_idx, train_idx] @ np.solve(M, y[:, train_idx])
            mse = np.power(ypred - y[:, test_idx], 2).mean(axis=0).T 
        '''
        ntrain = train_idx.shape[0]
        if subset_y_idx is None:
            subset_y_idx = np.arange(self.ny)
            
        M = (1 - theta_g) * np.eye(ntrain) + theta_g * self.grm[train_idx, train_idx]
        ypred = theta_g * (self.grm[test_idx, :][:, train_idx] @ np.linalg.solve(M, self.y[train_idx, :][:, subset_y_idx]))
        mse = np.power(ypred - self.y[test_idx, :][:, subset_y_idx], 2).mean(axis=0).T
        return mse, ypred
        
    def cv_train(self, train_idx, test_idx, rand_seed=1):
        np.random.seed(rand_seed)
        partitions = self.get_partition(train_idx)
        mse_mat = - np.ones((len(partitions), len(self.theta_g_grid), self.ny))
        for p in range(len(partitions)):
            inner_train_idx, inner_test_idx = partitions[p]
            for g in range(len(self.theta_g_grid)):
                theta_g = self.theta_g_grid[g]
                mse_mat[p, g, :], _ = self.train(inner_train_idx, inner_test_idx, theta_g)
        mse_mat = mse_mat.mean(axis=0) # squeeze partition dim by taking the average mse
        min_theta_idx = np.argmin(mse_mat, axis=0)
        y_pred = np.zeros((len(test_idx), self.ny))
        for g_cand in range(len(self.theta_g_grid)):
            theta_g = self.theta_g_grid[g_cand]
            y_active_idx = np.where(g_cand == min_theta_idx)[0]
            if y_active_idx.shape[0] == 0:
                continue
            _, y_active_pred = self.train(train_idx, test_idx, theta_g, subset_y_idx=y_active_idx)
            y_pred[:, y_active_idx] = y_active_pred
        return y_pred, self.y[test_idx, :]
            
