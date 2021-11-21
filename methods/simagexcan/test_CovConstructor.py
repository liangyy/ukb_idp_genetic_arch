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
def evd2mat(dict_):
    val = dict_['eig_val']
    vec = dict_['eig_vec']
    return vec @ np.diag(val) @ vec.T

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
constructor.compute_to_disk(
    mode='evd',
    param=0,
    output_prefix=output_prefix
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
res4 = np.load(f'{output_prefix}.evd.npz')
covmat4 = CovMatrix(f'{output_prefix}.evd.npz')
res4 = evd2mat(res4)

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
    
    # evd
    tmp = cov.copy()
    print('evd', np.allclose(res4, tmp))
    

print('---- testing cov eval_matmul_on_left ----')
x = np.random.rand(200, 19)
for cov in [ cov1, cov2 ]:
    
    # naive
    tmp = cov.copy()
    tmp = tmp @ x
    mul0, _ = covmat0.eval_matmul_on_left(x, param=2)
    print('naive', np.allclose(mul0, tmp))
    
    # cap
    tmp = cov.copy()
    tmp[ np.absolute(tmp) < threshold ] = 0
    tmp = tmp @ x
    mul1, _ = covmat1.eval_matmul_on_left(x)
    print('cap', np.allclose(mul1, tmp))
    
    # band 1
    tmp = get_band(cov.copy(), band1, tri=False)
    tmp = tmp @ x
    mul2, _ = covmat2.eval_matmul_on_left(x)
    print('band1', np.allclose(mul2, tmp))
    
    # band 2
    tmp = get_band(cov.copy(), band2, tri=False)
    tmp = tmp @ x
    mul3, _ = covmat3.eval_matmul_on_left(x)
    print('band2', np.allclose(mul3, tmp))
    
    # evd
    tmp = cov.copy()
    tmp = tmp @ x
    mul4, _ = covmat4.eval_matmul_on_left(x)
    print('evd', np.allclose(mul4, tmp))


print('---- testing cov eval_trace ----')
for cov in [ cov1, cov2 ]:
    
    # naive
    tmp = cov.copy()
    tmp = tmp.diagonal().sum()
    tr0 = covmat0.eval_trace(param=2)
    print('naive', np.allclose(tr0, tmp))
    
    # cap
    tmp = cov.copy()
    tmp[ np.absolute(tmp) < threshold ] = 0
    tmp = tmp.diagonal().sum()
    tr1 = covmat1.eval_trace()
    print('cap', np.allclose(tr1, tmp))
    
    # band 1
    tmp = get_band(cov.copy(), band1, tri=False)
    tmp = tmp.diagonal().sum()
    tr2 = covmat2.eval_trace()
    print('band1', np.allclose(tr2, tmp))
    
    # band 2
    tmp = get_band(cov.copy(), band2, tri=False)
    tmp = tmp.diagonal().sum()
    tr3 = covmat3.eval_trace()
    print('band2', np.allclose(tr3, tmp))
    
    # evd
    tmp = cov.copy()
    tmp = tmp.diagonal().sum()
    tr4 = covmat4.eval_trace()
    print('evd', np.allclose(tr4, tmp))    

print('---- testing cov eval_trace with cor=True ----')
for cov in [ cov1, cov2 ]:
    
    # naive
    tmp = cov.copy()
    d = np.sqrt(tmp.diagonal())
    tmp = tmp / d[:, np.newaxis] / d[np.newaxis, :]
    tmp = tmp.diagonal().sum()
    tr0 = covmat0.eval_trace(param=2, cor=True)
    print('naive', np.allclose(tr0, tmp))
    
    # cap
    tmp = cov.copy()
    tmp[ np.absolute(tmp) < threshold ] = 0
    d = np.sqrt(tmp.diagonal())
    tmp = tmp / d[:, np.newaxis] / d[np.newaxis, :]
    tmp = tmp.diagonal().sum()
    tr1 = covmat1.eval_trace(cor=True)
    print('cap', np.allclose(tr1, tmp))
    
    # band 1
    tmp = get_band(cov.copy(), band1, tri=False)
    d = np.sqrt(tmp.diagonal())
    tmp = tmp / d[:, np.newaxis] / d[np.newaxis, :]
    tmp = tmp.diagonal().sum()
    tr2 = covmat2.eval_trace(cor=True)
    print('band1', np.allclose(tr2, tmp))
    
    # band 2
    tmp = get_band(cov.copy(), band2, tri=False)
    d = np.sqrt(tmp.diagonal())
    tmp = tmp / d[:, np.newaxis] / d[np.newaxis, :]
    tmp = tmp.diagonal().sum()
    tr3 = covmat3.eval_trace(cor=True)
    print('band2', np.allclose(tr3, tmp))
    
    # not implemented
    # # evd
    # tmp = cov.copy()
    # tmp = tmp.diagonal().sum()
    # tr4 = covmat4.eval_trace()
    # print('evd', np.allclose(tr4, tmp))  


print('---- testing cov eval_sum_of_squares ----')
for cov in [ cov1, cov2 ]:
    
    # naive
    tmp = cov.copy()
    tmp = np.power(tmp, 2).sum()
    sq0 = covmat0.eval_sum_of_squares(param=2)
    print('naive', np.allclose(sq0, tmp))
    
    # cap
    tmp = cov.copy()
    tmp[ np.absolute(tmp) < threshold ] = 0
    tmp = np.power(tmp, 2).sum()
    sq1 = covmat1.eval_sum_of_squares()
    print('cap', np.allclose(sq1, tmp))
    
    # band 1
    tmp = get_band(cov.copy(), band1, tri=False)
    tmp = np.power(tmp, 2).sum()
    sq2 = covmat2.eval_sum_of_squares()
    print('band1', np.allclose(sq2, tmp))
    
    # band 2
    tmp = get_band(cov.copy(), band2, tri=False)
    tmp = np.power(tmp, 2).sum()
    sq3 = covmat3.eval_sum_of_squares()
    print('band2', np.allclose(sq3, tmp))
    
    # not implemented
    # # evd
    # tmp = cov.copy()
    # tmp = np.power(tmp, 2).sum()
    # tr4 = covmat4.eval_trace()
    # print('evd', np.allclose(tr4, tmp))   

print('---- testing cov eval_sum_of_squares with cor=True ----')
for cov in [ cov1, cov2 ]:
    
    # naive
    tmp = cov.copy()
    d = np.sqrt(tmp.diagonal())
    tmp = tmp / d[:, np.newaxis] / d[np.newaxis, :]
    tmp = np.power(tmp, 2).sum()
    sq0 = covmat0.eval_sum_of_squares(param=2, cor=True)
    print('naive', np.allclose(sq0, tmp))
    
    # cap
    tmp = cov.copy()
    tmp[ np.absolute(tmp) < threshold ] = 0
    d = np.sqrt(tmp.diagonal())
    tmp = tmp / d[:, np.newaxis] / d[np.newaxis, :]
    tmp = np.power(tmp, 2).sum()
    sq1 = covmat1.eval_sum_of_squares(cor=True)
    print('cap', np.allclose(sq1, tmp))
    
    # band 1
    tmp = get_band(cov.copy(), band1, tri=False)
    d = np.sqrt(tmp.diagonal())
    tmp = tmp / d[:, np.newaxis] / d[np.newaxis, :]
    tmp = np.power(tmp, 2).sum()
    sq2 = covmat2.eval_sum_of_squares(cor=True)
    print('band1', np.allclose(sq2, tmp))
    
    # band 2
    tmp = get_band(cov.copy(), band2, tri=False)
    d = np.sqrt(tmp.diagonal())
    tmp = tmp / d[:, np.newaxis] / d[np.newaxis, :]
    tmp = np.power(tmp, 2).sum()
    sq3 = covmat3.eval_sum_of_squares(cor=True)
    print('band2', np.allclose(sq3, tmp))
