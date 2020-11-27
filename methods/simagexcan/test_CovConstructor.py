import h5py
import numpy as np
from scipy.sparse import load_npz 
from CovConstructor import CovConstructor, CovMatrix

def get_band(mat, band_size, tri=True):
    i1 = np.tri(mat.shape[0], mat.shape[1], k=band_size)
    i2 = np.tri(mat.shape[0], mat.shape[1], k=-band_size - 1)
    mat[ np.logical_or(i1 == 0, i2 == 1) ] = 0
    if tri is True:
        return np.triu(mat)
    else:
        return mat

output_prefix = 'test_CovConstructor'

np.random.seed(2020)
mat = np.random.rand(100, 200)
threshold = 1e-3
band1 = 5
band2 = 40
constructor = CovConstructor(
    data=mat,
    nbatch=8
)
constructor.compute_to_disk(
    mode='naive',
    param=None,
    output_prefix=output_prefix
)
constructor.compute_to_disk(
    mode='cap',
    param=threshold,
    output_prefix=output_prefix
)
constructor.compute_to_disk(
    mode='banded',
    param=band1,
    output_prefix=output_prefix
)
constructor.compute_to_disk(
    mode='banded',
    param=band2,
    output_prefix=output_prefix + '2'
)
with h5py.File(f'{output_prefix}.naive.h5', 'r') as f:
    res0 = f['cov'][:]
covmat0 = CovMatrix(f'{output_prefix}.naive.h5')
res1 = load_npz(f'{output_prefix}.cap.npz').todense()
covmat1 = CovMatrix(f'{output_prefix}.cap.npz')
res2 = load_npz(f'{output_prefix}.banded.npz').todense()
covmat2 = CovMatrix(f'{output_prefix}.banded.npz')
res3 = load_npz(f'{output_prefix}2.banded.npz').todense()
covmat3 = CovMatrix(f'{output_prefix}2.banded.npz')

cov1 = np.cov(mat.T)
mat_centered = mat - mat.mean(axis=0)
cov2 = mat_centered.T @ mat_centered / (mat_centered.shape[0] - 1)

print('---- testing cov construction ----')
for cov in [ cov1, cov2 ]:
    
    # naive
    tmp = cov.copy()
    tmp = np.triu(tmp)
    print('naive', np.allclose(res0, tmp))
    
    # cap
    tmp = cov.copy()
    tmp = np.triu(tmp)
    tmp[ np.absolute(tmp) < threshold ] = 0
    print('cap', np.allclose(res1, tmp))
    
    # band 1
    tmp = get_band(cov.copy(), band1)
    print('band1', np.allclose(res2, tmp))
    
    # band 2
    tmp = get_band(cov.copy(), band2)
    print('band2', np.allclose(res3, tmp))
    

print('---- testing cov eval_matmul_on_left ----')
x = np.random.rand(200, 19)
for cov in [ cov1, cov2 ]:
    
    # naive
    tmp = cov.copy()
    tmp = tmp @ x
    mul0 = covmat0.eval_matmul_on_left(x)
    print('naive', np.allclose(mul0, tmp))
    
    # cap
    tmp = cov.copy()
    tmp[ np.absolute(tmp) < threshold ] = 0
    tmp = tmp @ x
    mul1 = covmat1.eval_matmul_on_left(x)
    print('cap', np.allclose(mul1, tmp))
    
    # band 1
    tmp = get_band(cov.copy(), band1, tri=False)
    tmp = tmp @ x
    mul2 = covmat2.eval_matmul_on_left(x)
    print('band1', np.allclose(mul2, tmp))
    
    # band 2
    tmp = get_band(cov.copy(), band2, tri=False)
    tmp = tmp @ x
    mul3 = covmat3.eval_matmul_on_left(x)
    print('band2', np.allclose(mul3, tmp))
    