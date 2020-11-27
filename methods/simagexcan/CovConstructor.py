import numpy as np
from scipy.sparse import coo_matrix, save_npz

class CovConstructor:
    def __init__(self, data, nbatch=None):
        '''
        Important: self.data is ALWAYS centered
        '''
        self.data = data - data.mean(axis=0)
        self.ncol = self.data.shape[1]
        self.nrow = self.data.shape[0]
        self.nbatch = nbatch
        self._set_batch_anchors()
    def _set_batch_anchors(self):
        ncol = self.ncol
        batch_size = ncol // self.nbatch
        if batch_size * self.nbatch < ncol:
            batch_size += 1
        if batch_size < 10:
            raise ValueError('Too many batches. Exit.')
        self.batch_anchors = []
        for i in range(self.nbatch):
            start = i * batch_size
            end = min((i + 1) * batch_size, ncol)
            self.batch_anchors.append([start, end])
    def _flatten_2d_mat(self, mat):
        row_index_mat = np.tile(np.arange(mat.shape[0]), reps=(mat.shape[1], 1)).T
        row = row_index_mat.flatten()
        del row_index_mat
        col_index_mat = np.tile(np.arange(mat.shape[1]), reps=(mat.shape[0], 1))
        col = col_index_mat.flatten()
        del col_index_mat
        return row, col, mat.flatten()
    def compute_to_h5(fn):
        raise NotImplementedError('Not yet implemented `compute_to_h5`.')
    def compute_to_disk(self, mode, output_prefix, param=None):
        if mode == 'naive':
            fn = output_prefix + '.h5'
            self.compute_to_h5(fn)
        elif mode == 'cap':
            fn = output_prefix + '.cap.npz'
            self.compute_to_cap_npz(fn, threshold=param)
        elif mode == 'banded':
            fn = output_prefix + '.banded.npz'
            self.compute_to_banded_npz(fn, band_size=param)
    def _compute_cov(self, s1, e1, s2, e2):
        '''
        Given submatrix index: 
            matrix1 = [:, s1 : e1], matrix2 = [:, s2 : e2]
        Return: 
            Pairwise covariance between column in matrix1 and column in matrix2.
            Elements are returned in row, col, val format. 
            And only row <= col ones are returned.
        Formula:
            covariance = col1 * col2 / self.nrow (col is centered in __init__) 
        '''
        tmp = np.einsum('ni,nj->ij', self.data[:, s1 : e1], self.data[:, s2 : e2]) / (self.nrow - 1)
        row, col, val = self._flatten_2d_mat(tmp)
        row += s1
        col += s2
        to_keep = row <= col
        row, col, val = row[to_keep], col[to_keep], val[to_keep]
        del tmp
        return row, col, val
        
    def compute_to_banded_npz(self, fn, band_size=100):
        row_all, col_all, value_all = [], [], []
        for i, (s1, e1) in enumerate(self.batch_anchors):
            for j, (s2, e2) in enumerate(self.batch_anchors):
                if i > j:
                    continue
                if s2 > e1 + band_size - 1:
                    continue
                row, col, value = self._compute_cov(s1, e1, s2, e2)
                to_keep = col - row <= band_size 
                row, col, value = row[to_keep], col[to_keep], value[to_keep]
                row_all.append(row)
                col_all.append(col)
                value_all.append(value)
        row_all = np.concatenate(row_all, axis=0)
        col_all = np.concatenate(col_all, axis=0)
        value_all = np.concatenate(value_all, axis=0)
        cov_coo = coo_matrix(
            (value_all, (row_all, col_all)), 
            shape=(self.ncol, self.ncol)
        )
        save_npz(fn, cov_coo)  
    def compute_to_cap_npz(self, fn, threshold=1e-5):
        row_all, col_all, value_all = [], [], []
        for i, (s1, e1) in enumerate(self.batch_anchors):
            for j, (s2, e2) in enumerate(self.batch_anchors):
                if i > j:
                    continue
                row, col, value = self._compute_cov(s1, e1, s2, e2)
                to_keep = np.abs(value) > threshold
                row, col, value = row[to_keep], col[to_keep], value[to_keep]
                row_all.append(row)
                col_all.append(col)
                value_all.append(value)
        row_all = np.concatenate(row_all, axis=0)
        col_all = np.concatenate(col_all, axis=0)
        value_all = np.concatenate(value_all, axis=0)
        cov_coo = coo_matrix(
            (value_all, (row_all, col_all)), 
            shape=(self.ncol, self.ncol)
        )
        save_npz(fn, cov_coo)            
            
            