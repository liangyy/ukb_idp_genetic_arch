import numpy as np
from scipy.sparse import coo_matrix, save_npz, load_npz, csr_matrix

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
        if batch_size < 10 and self.nbatch > 5:
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
    def compute_to_h5(self, fn, dtype='f'):
        import h5py
        f = h5py.File(fn, 'w')
        dset = f.create_dataset('cov', (self.ncol, self.ncol), dtype=dtype)
        for i, (s1, e1) in enumerate(self.batch_anchors):
            for j, (s2, e2) in enumerate(self.batch_anchors):
                if i > j:
                    continue
                dset[s1 : e1, s2 : e2] = self._compute_cov(
                    s1, e1, s2, e2, 
                    flatten=False, triu=(i == j)
                )
                tmp = self._compute_cov(
                    s1, e1, s2, e2, 
                    flatten=False, triu=(i == j)
                )
        f.close()
    def compute_to_disk(self, mode, output_prefix, param=None):
        if mode == 'naive':
            fn = output_prefix + '.naive.h5'
            self.compute_to_h5(fn, dtype=param)
        elif mode == 'cap':
            fn = output_prefix + '.cap.npz'
            self.compute_to_cap_npz(fn, threshold=param)
        elif mode == 'banded':
            fn = output_prefix + '.banded.npz'
            self.compute_to_banded_npz(fn, band_size=param)
        elif mode == 'evd':
            fn = output_prefix + '.evd.npz'
            self.compute_to_evd_npz(fn, min_max_thres=param)
    def _compute_cov(self, s1, e1, s2, e2, flatten=True, triu=True):
        '''
        Given submatrix index: 
            matrix1 = [:, s1 : e1], matrix2 = [:, s2 : e2]
        Return: 
            Pairwise covariance between column in matrix1 and column in matrix2.
            Elements are returned in row, col, val format (flatten = True). 
            And only row <= col ones are returned.
            But if flatten = False, triu could be set to False
            to return the full matrix.
        Formula:
            covariance = col1 * col2 / self.nrow (col is centered in __init__) 
        '''
        tmp = np.einsum('ni,nj->ij', self.data[:, s1 : e1], self.data[:, s2 : e2]) / (self.nrow - 1)
        if flatten is False:
            if triu is True:
                return np.triu(tmp)
            else:
                return tmp
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
    def compute_to_evd_npz(self, fn, min_max_thres=0):
        
        # determine mode
        if self.ncol > self.nrow:
            mode = 'XXt'
        else:
            mode = 'XtX'
        
        # run evd
        if mode == 'XtX':
            target_mat = self.data.T @ self.data / (self.nrow - 1)
        elif mode == 'XXt':
            target_mat = self.data @ self.data.T / (self.nrow - 1)
        
        eig_val, eig_vec = np.linalg.eigh(target_mat)
        
        # thresholding
        max_ = eig_val[-1]
        if max_ < 0:
            raise ValueError('Largest eigen value is smaller than zero. Something wrong.')
        to_keep_ind = eig_val / max_ > min_max_thres
        eig_vec = eig_vec[:, to_keep_ind]
        eig_val = eig_val[to_keep_ind]
        
        # calculate V for XXt
        if mode == 'XXt':
            # XX' / (n - 1) = U S U'
            # X'X / (n - 1) = V S V'
            # L^-1 U' X / sqrt(n - 1) = L^-1 U' U L V' = V'
            # L = sqrt(S)
            eig_vec = self.data.T @ eig_vec / np.sqrt(eig_val)[np.newaxis, :] / np.sqrt(self.nrow - 1)
        
        # save to disk
        np.savez(fn, eig_val=eig_val, eig_vec=eig_vec)  
        
class CovMatrix:
    def __init__(self, fn):
        self.fn = fn
        self.mode = self._init_mode()
    def _init_mode(self):
        import pathlib
        if not pathlib.Path(self.fn).is_file():
            raise ValueError('Input file does not exist.') 
        tmp = self.fn.split('.')
        return tmp[-2]
    def eval_matmul_on_left(self, left_mat, param=None):
        '''
        Retur cov @ left_mat along with the diag 
        '''
        if self.mode in ['banded', 'cap']:
            return self._eval_matmul_on_left_npz(left_mat)
        elif self.mode == 'naive':
            return self._eval_matmul_on_left_h5(left_mat, batch_size=param)
        elif self.mode == 'evd':
            return self._eval_matmul_on_left_evd(left_mat)
    def eval_trace(self, param=None, cor=False):
        if self.mode in ['banded', 'cap']:
            return self._eval_trace_npz(cor=cor)
        elif self.mode == 'naive':
            return self._eval_trace_h5(batch_size=param, cor=cor)
        elif self.mode == 'evd':
            return self._eval_trace_evd(cor=cor)
    def eval_sum_of_squares(self, param=None, cor=False):
        '''
        Evaluate sum_{ij} R_{ij}^2
        '''
        if self.mode in ['banded', 'cap']:
            return self._eval_sum_of_squares_npz(cor=cor)
        elif self.mode == 'naive':
            return self._eval_sum_of_squares_h5(batch_size=param, cor=cor)
        elif self.mode == 'evd':
            return self._eval_sum_of_squares_evd(cor=cor)
    @staticmethod
    def _npz_standardize(csr, cor):
        if cor is True:
            diag_sqrt = np.sqrt(csr.diagonal())
            # adapted from https://stackoverflow.com/questions/49254111/row-division-in-scipy-sparse-matrix
            r, c = csr.nonzero()
            row_sp = csr_matrix(
                ((1.0 / diag_sqrt)[r], (r,c)), 
                shape=(csr.shape))
            col_sp = csr_matrix(
                ((1.0 / diag_sqrt)[c], (r,c)), 
                shape=(csr.shape))
            csr = csr.multiply(row_sp)
            csr = csr.multiply(col_sp)
        return csr
    def _eval_sum_of_squares_evd(self, cor):
        raise NotImplementedError(
            'Have not implemented sum of squres for evd cov')
    def _eval_sum_of_squares_h5(self, batch_size, cor):
        import h5py
        f = h5py.File(self.fn, 'r')
        nrow = f['cov'].shape[0]
        if batch_size is None:
            nbatch = 1
            batch_size = nrow
        else:
            nbatch = nrow // batch_size
            if nbatch * batch_size < nrow:
                nbatch += 1
        s, e = 0, batch_size
        upper = np.zeros(nrow)
        diag_vec = np.zeros(nrow)
        for i in range(nbatch):
            tmp = f['cov'][s : e, s : ]
            diag_i = tmp.diagonal()
            diag_vec[s : e] = diag_i
            if cor is True:
                tmp = tmp / np.sqrt(diag_i[:, np.newaxis])
            upper[s : ] += np.power(np.triu(tmp, 1), 2).sum(axis=0) 
            s += batch_size
            e = min(batch_size + e, nrow)
        f.close()
        if cor is False:
            diag = np.power(diag_vec, 2).sum()
        else:
            upper = upper / diag_vec
            diag = nrow
        return upper.sum() * 2 + diag
    def _eval_sum_of_squares_npz(self, cor):
        csr = load_npz(self.fn).tocsr()
        csr = self._npz_standardize(csr, cor)
        return csr.power(2).sum() * 2 - np.power(csr.diagonal(), 2).sum()
    def _eval_trace_npz(self, cor):
        csr = load_npz(self.fn).tocsr()
        csr = self._npz_standardize(csr, cor)
        return csr.diagonal().sum()
    def _eval_trace_h5(self, batch_size, cor):
        import h5py
        f = h5py.File(self.fn, 'r')
        nrow = f['cov'].shape[0]
        if cor is True:
            return nrow
        if batch_size is None:
            nbatch = 1
            batch_size = nrow
        else:
            nbatch = nrow // batch_size
            if nbatch * batch_size < nrow:
                nbatch += 1
        s, e = 0, batch_size
        diag_cov = np.zeros((nrow))
        for i in range(nbatch):
            diag_cov[s : e] = f['cov'][s : e, s : e].diagonal()
            s += batch_size
            e = min(batch_size + e, nrow)
        f.close()
        return diag_cov.sum()
    def _eval_trace_evd(self, cor):
        if cor is True:
            raise NotImplementedError(
                'Have not implemented trace for evd cor')
        res = np.load(self.fn)
        if 'eig_val' not in res or 'eig_vec' not in res:
            raise ValueError('There are {} in file but missing either eig_val or eig_vec.'.format(list(res.keys())))
        eig_val, eig_vec = res['eig_val'], res['eig_vec']
        return np.einsum('ij,j,ij', eig_vec, eig_val, eig_vec)
    def _eval_matmul_on_left_evd(self, mat):
        res = np.load(self.fn)
        if 'eig_val' not in res or 'eig_vec' not in res:
            raise ValueError('There are {} in file but missing either eig_val or eig_vec.'.format(list(res.keys())))
        eig_val, eig_vec = res['eig_val'], res['eig_vec']
        s1 = eig_vec.T @ mat
        s2 = eig_val[:, np.newaxis] * s1
        s3 = eig_vec @ s2
        s4 = np.einsum('ij,j,ij->i', eig_vec, eig_val, eig_vec)
        return s3, s4 
    def _eval_matmul_on_left_npz(self, mat):
        csr = load_npz(self.fn).tocsr()
        diag_csr = csr.diagonal()
        return csr.dot(mat) + csr.transpose().dot(mat) - diag_csr[:, np.newaxis] * mat, diag_csr
    def _eval_matmul_on_left_h5(self, mat, batch_size=None):
        import h5py
        f = h5py.File(self.fn, 'r')
        nrow = f['cov'].shape[0]
        ncol = mat.shape[1]
        res = np.zeros((nrow, ncol))
        if batch_size is None:
            nbatch = 1
            batch_size = nrow
        else:
            nbatch = nrow // batch_size
            if nbatch * batch_size < nrow:
                nbatch += 1
        s, e = 0, batch_size
        diag_cov = np.zeros((nrow))
        for i in range(nbatch):
            res[s : e, :] = f['cov'][s : e, :] @ mat
            res[s : e, :] += f['cov'][:, s : e].T @ mat
            diag_cov[s : e] = f['cov'][s : e, s : e].diagonal()
            s += batch_size
            e = min(batch_size + e, nrow)
        res -= diag_cov[:, np.newaxis] * mat
        f.close()
        return res, diag_cov
        
