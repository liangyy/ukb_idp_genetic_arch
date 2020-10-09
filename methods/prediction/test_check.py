import pandas as pd
import numpy as np
e2 = pd.read_csv('test_prs.sscore', sep='\t')
e = pd.read_parquet('test_prs.parquet')
e.indiv = e.indiv.astype(int); e = pd.merge(e, e2, left_on='indiv', right_on='#IID')
ll = [ 'IDP-25303', 'IDP-25304', 'IDP-25305' ]
for xx in ll:
    for yy in ll:
        print(xx, 'and', yy, np.corrcoef(e[xx], e[yy + '_SUM'])[0,1])
# print(np.corrcoef(e['IDP-25304'], e['IDP-25304_SUM']))
# print(np.corrcoef(e['IDP-25305'], e['IDP-25305_SUM']))

