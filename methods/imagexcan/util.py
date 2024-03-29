import yaml
import os
import numpy as np
import pandas as pd
import scipy.stats

# need to export PYTHONPATH to include https://github.com/liangyy/misc-tools/blob/master/pyutil/
from pystat import inv_norm_vec, standardize_by_col

def bhat_pval_to_zscore(bhat, pval):
    z_abs = scipy.stats.norm.isf(pval / 2)
    return z_abs * np.sign(bhat)

def calc_cor(x, covar=None):
    if covar is not None:
        tmpq, tmpr = np.linalg.qr(covar, mode='reduced')
        x_std = x - tmpq @ (tmpq.T @ x)
        x_std = standardize_by_col(x_std)
    else:
        x_std = standardize_by_col(x)
    return (x_std.T @ x_std) / x_std.shape[0]

def read_yaml(yaml_):
    with open(yaml_, 'r') as f:
        o = yaml.safe_load(f)
    return o

def load_list(fn):
    o = []
    with open(fn, 'r') as f:
        for i in f:
            o.append(i.strip())
    return o

def take_intersect(lists):
    if len(lists) == 0:
        return lists 
    o = set(lists[0])
    for l in lists:
        o = o.intersection(set(l))
    return list(o)

def exclude_b_from_a(a, b):
    a_ = set(a)
    b_ = set(b)
    a_ = a_.difference(b_)
    return list(a_)

def check_a_in_b(a, b):
    return all(i in b for i in a)

def check_binary(x):
    xr = np.round_(x)
    return np.logical_not(np.logical_or(xr == 1, xr == 0)).sum() == 0

def rearrange_rows(df, target_list):
    tmp = pd.DataFrame({'indiv': target_list})
    return pd.merge(tmp, df, on='indiv', how='left')

def get_matrix(df):
    '''
    Assume first col is individual id.
    Return the values without individual id.
    '''
    if df is None:
        return None
    else:
        return df.iloc[:, 1:].values

def load_table(file_list):
    fn, indiv_col = file_list[:2]
    _, fn_ext = os.path.splitext(fn)
    if fn_ext == '.parquet':
        df = pd.read_parquet(fn)
    elif fn_ext == '.csv':
        df = pd.read_csv(fn)
    elif fn_ext == '.txt' or fn_ext == '.tsv':
        df = pd.read_csv(fn, sep='\s+')
    for i in range(df.shape[1]):
        if df.columns[i] == indiv_col:
            break
    col_list = df.columns.to_list()
    col_list.pop(i)
    col_list = [ indiv_col ] + col_list
    df = df.reindex(columns=col_list)
    df.rename(columns={indiv_col: 'indiv'}, inplace=True)
    df.indiv = df.indiv.astype(str)
    return df

def load_idp(f_list):
    return load_table(f_list)

def load_phenotype(file_list, fyaml):
    df = load_table(file_list)
    df_yaml = read_yaml(fyaml)
    if not check_a_in_b(list(df_yaml.keys()), df.columns[1:].to_list()):
        raise ValueError(f'The {file_list} does not match {fyaml}.')
    df = df.reindex(columns=['indiv'] + list(df_yaml.keys()))
    df.drop_duplicates(subset='indiv', inplace=True)
    return df, df_yaml

def load_covariate(file_list, fyaml):
    df = load_table(file_list)
    df_yaml = read_yaml(fyaml)
    if not check_a_in_b(list(df_yaml.keys()), df.columns.to_list()):
        raise ValueError(f'The {file_string} does not match {fyaml}.')
    df_res = [ pd.DataFrame({'indiv': df.indiv}) ]
    nsample = df.shape[0]
    for col, ctype in df_yaml.items():
        if ctype == 'continuous':
            df_res.append(pd.DataFrame({col: df[col]}))
        if ctype == 'categorical':
            col_values = df[col].values
            categories = np.unique(col_values)
            ncate = categories.shape[0]
            if ncate > 25:
                raise ValueError(f'Too many categories at {col}.')
            mat = np.zeros((nsample, ncate))
            for i in range(ncate):
                mat[ col_values == categories[i], i ] = 1
            df_res.append(
                pd.DataFrame(
                    mat, 
                    columns=[ f'{col}_{val}' for val in categories.tolist() ]
                )
            )
    return pd.concat(df_res, axis=1)

def intersection(l1, l2):
    a = set(l1)
    a = a.intersection(set(l2))
    return sorted(list(a))
    
def to_list(var):
    if not isinstance(var, list):
        return [ var ]
    else:
        return var
