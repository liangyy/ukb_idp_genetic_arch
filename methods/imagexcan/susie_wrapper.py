import pandas as pd

import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import numpy2ri
numpy2ri.activate()
r_base = importr('base')
r_susieR = importr('susieR')

def run_susie_wrapper(z, R, params={}):# kwargs={}):
    mod = r_susieR.susie_rss(r_base.as_vector(z), R, **params)
    df_ = r_base.summary(mod).rx2('vars')
    with localconverter(ro.default_converter + pandas2ri.converter):
        df_ = ro.conversion.rpy2py(df_)
    df_ = pd.DataFrame({'pip': list(df_.variable_prob), 'cs': list(df_.cs), 'idx': list(df_.variable)})
    df_.sort_values(by='idx', inplace=True)
    return list(df_.pip), list(df_.cs)