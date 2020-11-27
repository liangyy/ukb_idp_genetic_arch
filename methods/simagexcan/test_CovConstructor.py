import numpy as np
from scipy.sparse import load_npz 
from CovConstructor import CovConstructor

def get_band(mat, band_size):
    i1 = np.tri(mat.shape[0], mat.shape[1], k=band_size)
    i2 = np.tri(mat.shape[0], mat.shape[1], k=-band_size)
    mat[ np.logical_or(i1 == 0, i2 == 1) ] = 0
    return np.triu(mat)

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
res1 = load_npz(f'{output_prefix}.cap.npz').todense()
res2 = load_npz(f'{output_prefix}.banded.npz').todense()
res3 = load_npz(f'{output_prefix}2.banded.npz').todense()
cov1 = np.cov(mat.T)
mat_centered = mat - mat.mean(axis=0)
cov2 = mat_centered.T @ mat_centered / (mat_centered.shape[0] - 1)


for cov in [ cov1, cov2 ]:
    
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
    
    