import sys
import pandas as pd
import numpy as np
import yaml
from collections import OrderedDict 
np.random.seed(2020)

file_in = sys.argv[1]
file_out1 = sys.argv[2]
file_out2 = sys.argv[3]
file_out3 = sys.argv[4]

K = 5

df = pd.read_parquet(file_in)
rand_idps = np.random.randint(0, df.shape[1] - 1, size=K) + 1  # the offset is indiv col
pheno_names = []
out = OrderedDict()
out['indiv'] = df.indiv
for i, k in enumerate(list(rand_idps)):
    name = '{}_{}'.format(i, df.columns[k])
    tmp = df.iloc[:, k].values
    tmp = tmp - tmp.mean()
    out[name] = tmp + np.random.normal(scale=4 * tmp.std(), size=df.shape[0]) + np.random.randint(-5, 5, size=1)
    name = '{}_null'.format(i, df.columns[k])
    out[name] = np.random.normal(scale=4 * tmp.std(), size=df.shape[0]) + np.random.randint(-5, 5, size=1)
out = pd.DataFrame(out)
out.to_csv(file_out1, index=False)

with open(file_out2, 'w') as f:
    yaml.dump({ k: 'linear_regression' for k in out.keys()[1:] }, f)

with open(file_out3, 'w') as f:
    yaml.dump({ k: 'susie' for k in out.keys()[1:] }, f)
    

