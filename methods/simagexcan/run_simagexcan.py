import re
import pandas as pd
import numpy as np
import pathlib
from collections import OrderedDict 
from pyutil import read_table

BASE_PAIR = {
    'A': 'T',
    'T': 'A',
    'G': 'C',
    'C': 'G'
}

def check_flip(a1, a2, b1, b2):
    res = []
    for _a1, _a2, _b1, _b2 in zip(a1, a2, b1, b2):
        res.append(_check_flip(_a1, _a2, _b1, _b2))
    return np.array(res)

def _check_flip(a0, a1, b0, b1):
    '''
    check if (a0, a1) and (b0, b1) are of the same direction.
    If there is nan or they don't match at all or ambiguious return nan
    Else if they are in the same direction, return 1
    Else return -1
    '''
    if a0 is np.nan or a1 is np.nan or b0 is np.nan or b1 is np.nan:
        return np.nan
    # remove ambiguious first.
    if a0 == BASE_PAIR[a1] or b0 == BASE_PAIR[b1]:
        return np.nan
    # exact match
    if a0 == b0 and a1 == b1:
        return 1
    # flip
    if a0 == b1 and a1 == b0:
        return -1    
    # compliment match
    if a0 == BASE_PAIR[b0] and a1 == BASE_PAIR[b1]:
        return 1
    # compliment flip
    if a0 == BASE_PAIR[b1] and a1 == BASE_PAIR[b0]:
        return -1  
    # if all above does not return, it has to be invalid.
    return np.nan

def rearrage_df_by_target(df, target, df_value_cols):
    df_res = target[['snpid', 'chr', 'effect_allele', 'non_effect_allele']]
    df_res = pd.merge(
        df_res, df, 
        on=['snpid', 'chr'], 
        suffixes=['_res', '_df'],
        how='left'
    )
    flip_factor = check_flip(
        a1=df_res.effect_allele_res, 
        a2=df_res.non_effect_allele_res,
        b1=df_res.effect_allele_df, 
        b2=df_res.non_effect_allele_df
    )
    # we need to carry the missingness when we move on
    with np.errstate(invalid='ignore'):
        df_res[df_value_cols] = df_res[df_value_cols] * flip_factor[:, np.newaxis]
    df_res.drop(
        columns=['effect_allele_df', 'non_effect_allele_df'], inplace=True
    )
    df_res.rename(
        columns={
            'effect_allele_res': 'effect_allele',
            'non_effect_allele_res': 'non_effect_allele'
        }, 
        inplace=True
    )
    return df_res

def harmonize_gwas_and_weight(gwas, weight):
    '''
    Harmonize GWAS to weight SNP set.
    But only keep the ones that present in both.
    '''
    df_common = pd.merge(
        gwas[['snpid', 'chr', 'effect_allele', 'non_effect_allele']],
        weight[['snpid', 'chr', 'effect_allele', 'non_effect_allele']],
        on=['snpid', 'chr'],
        suffixes=['_gwas', '_weight']
    )
    flip_factor = check_flip(
        a1=df_common.effect_allele_gwas, 
        a2=df_common.non_effect_allele_gwas,
        b1=df_common.effect_allele_weight, 
        b2=df_common.non_effect_allele_weight
    )
    
    # need to remove the invalid variant before moving on
    to_keep_ind = np.logical_not(np.isnan(flip_factor))
    df_common = df_common[ to_keep_ind ].reset_index(drop=True)
    flip_factor = flip_factor[ to_keep_ind ]
    
    df_common.drop(columns=['effect_allele_gwas', 'non_effect_allele_gwas'], inplace=True)
    df_common.rename(columns={'effect_allele_weight': 'effect_allele', 'non_effect_allele_weight': 'non_effect_allele'}, inplace=True)

    df_gwas = pd.merge(
        df_common[['snpid', 'chr', 'effect_allele', 'non_effect_allele']], 
        gwas.drop(columns=['effect_allele', 'non_effect_allele']), 
        on=['snpid', 'chr']
    )
    df_gwas.effect_size = df_gwas.effect_size * flip_factor
    
    df_weight = pd.merge(
        df_common[['snpid', 'chr', 'effect_allele', 'non_effect_allele']], 
        weight.drop(columns=['effect_allele', 'non_effect_allele']), 
        on=['snpid', 'chr']
    )
    return df_gwas, df_weight

def _parse_args(args_list, desired_cols):
    fn = args_list[0]
    if not pathlib.Path(fn).is_file():
        raise ValueError('Filename is wrong. Cannot find the file.')
    dict = {}
    snpid_name = None
    for i in args_list[1:]:
        tmp = i.split(':')
        if len(tmp) != 2:
            raise ValueError('Wrong gwas args list. Need [col]:[name] pairs.')
        col, name = tmp
        if col not in desired_cols:
            raise ValueError(f'Wrong col = {col}.')
        dict[col] = name
    rename_dict = OrderedDict()
    for dd in desired_cols:
        if dd not in dict:
            raise ValueError(f'Need to have col = {dd}.')
        rename_dict[dict[dd]] = dd
    return fn, rename_dict
    
def _parse_gwas_args(args_list):
    have_effect_size = False
    for kk in args_list:
        if 'effect_size:' in kk:
            have_effect_size = True
    if have_effect_size is True:
        desired_cols = [
            'snpid', 'non_effect_allele', 'effect_allele', 
            'effect_size', 'effect_size_se', 'chr'
        ]
    else:
        desired_cols = [
            'snpid', 'non_effect_allele', 'effect_allele', 
            'zscore', 'allele_frequency', 'sample_size', 'chr'
        ]
    fn, rename_dict = _parse_args(args_list, desired_cols)
    for k, v in rename_dict.items():
        if v == 'snpid':
            snpid_name = k
            break
    return fn, rename_dict, snpid_name

def impute_b_from_z(zscore, af, n):
    se = 1 / np.sqrt(2 * n * af * (1 - af))
    bhat = zscore * se
    return bhat, se

def clean_up_chr(ll):
    for i in range(len(ll)):
        ll[i] = re.sub('chr', '', ll[i])
    return ll
    
def load_gwas(gwas_args_list):
    fn, rename_dict, snpid_col = _parse_gwas_args(gwas_args_list)
    df = read_table(fn, indiv_col=snpid_col)
    df.rename(columns={'indiv': snpid_col}, inplace=True)
    df.rename(columns=rename_dict, inplace=True)
    df.drop_duplicates('snpid', inplace=True)
    df.chr = clean_up_chr(list(df.chr.astype(str)))
    if 'effect_size' not in rename_dict.values():
        df['effect_size'], df['effect_size_se'] = impute_b_from_z(df.zscore, df.allele_frequency, df.sample_size)
    desired_cols = [
        'snpid', 'non_effect_allele', 'effect_allele', 
        'effect_size', 'effect_size_se', 'chr'
    ]
    return df[desired_cols]

def _parse_idp_args(args_list):
    desired_cols = [
        'snpid', 'non_effect_allele', 'effect_allele', 'chr'
    ]
    fn, rename_dict = _parse_args(args_list, desired_cols)
    return fn, rename_dict
        
def load_idp(args_list):
    fn, rename_dict = _parse_idp_args(args_list)
    df = pd.read_parquet(fn)
    df.rename(columns=rename_dict, inplace=True)
    df.chr = df.chr.astype(str)
    return df

def load_cov_meta(fn):
    fn = '.'.join(fn.split('.')[:-1])
    fn = fn + '.snp_meta.parquet'
    return pd.read_parquet(fn)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='run_simagexcan.py', description='''
        Run S-ImageXcan with pre-computed genotype covariance.
        Need to export PYTHONPATH=path-to/imagexcan:path-to/misc-tools/pyutil
    ''')
    parser.add_argument('--genotype_covariance', help='''
        The genotype covariance computed in build_genotype_covariance.py
        Accept wildcard {chr_num}.
        Will automatically search for the corresponding meta SNP file.
    ''')
    parser.add_argument('--gwas', nargs='+', help='''
        Need to have column names for: 
            snpid, non_effect_allele, effect_allele, 
            effect_size, effect_size_se, chr.
        If there is no effect_size avaliable, it could 
        impute effect_size from zscore, allele_frequency, 
        sample_size.
        The format is: snpid:rsid_col, ..., chr:chr
    ''')
    parser.add_argument('--idp_weight', nargs='+', help='''
        The IDP weight table is in parquet format.
        It contains columns:
            snpid, effect_allele, non_effect_allele, chr.
        Along with all other columns for the IDPs.
        Specify the column names, e.g.: snpid:rsID, ..., chr:chr
    ''')
    parser.add_argument('--output', help='''
        The output CSV filename.
        Will return both marginal test result and also the susieR result.
    ''')
    parser.add_argument('--z_ld_weight', type=float, default=1e-4, help='''
        LD = (1 - z_ld_weight) * LD + z_ld_weight * (Z @ Z.T)
        to avoid mis-specified LD.
    ''')
    args = parser.parse_args()
    
    from tqdm import tqdm
    import logging, time, sys, os
    # configing util
    logging.basicConfig(
        level = logging.INFO, 
        stream = sys.stderr, 
        format = '%(asctime)s  %(message)s',
        datefmt = '%Y-%m-%d %I:%M:%S %p'
    )
    from CovConstructor import CovMatrix
    from susie_wrapper import run_susie_wrapper
    from pystat import z2p
    
    logging.info('Loading GWAS.')
    df_gwas = load_gwas(args.gwas)
    # df_gwas columns: 
    # snpid, non_effect_allele, effect_allele, 
    # effect_size, effect_size_se, chr
    logging.info('GWAS SNP = {}'.format(df_gwas.shape[0]))
    
    logging.info('Loading IDP weights.')
    df_weight = load_idp(args.idp_weight)
    idp_names = list(df_weight.columns[4:])
    nidp = len(idp_names)
    logging.info('IDP SNP = {} and number of IDPs = {}'.format(df_weight.shape[0], nidp))
    
    logging.info('Harmonizing GWAS and IDP weights.')
    # harmonize GWAS and IDP weight table so that they have the same set of 
    # SNPs (including direction).
    df_gwas, df_weight = harmonize_gwas_and_weight(df_gwas, df_weight)
    logging.info('{} SNPs left after harmonizing GWAS and IDP weights.'.format(df_gwas.shape[0]))
    
    # please refer to https://github.com/hakyimlab/yanyu-notebook/blob/master/notes/date_112420.Rmd
    # for the details of the S-ImageXcan formula
    # to take the following procedure.
    # 0. subset IDP and GWAS SNPs.
    # 1. Per chromosome
    #   1.1 obtain D(chr), S_R(chr), and var_R(chr).
    #   1.2 compute numer_b(chr) = Gamma(chr).T @ (var_R(chr) * b_gwas(chr)) 
    #   1.3 compute numer_z(chr) = Gamma(chr).T @ (S_R(chr) * z_gwas(chr)) 
    # 2. compute marginal test.
    #   2.1 D = sum_chr D(chr), var_D = diag(D), S_D = sqrt(var_D)
    #   2.2 beta_imagexcan = ( sum_chr numer_b(chr) ) / var_D
    #   2.3 z_imagexcan =  ( sum_chr numer_z(chr) ) / S_D
    # 3. run susieR.
    #   3.1 Sigma = D / S_D[:, np.newaxis] / S_D[np.newaxis, :]
    
    D = np.zeros((nidp, nidp))
    numer_b = np.zeros((nidp))
    numer_z = np.zeros((nidp))
    for i in range(1, 23):
        
        df_gwas_sub = df_gwas[ df_gwas.chr == str(i) ].reset_index(drop=True)
        df_weight_sub = df_weight[ df_weight.chr == str(i) ].reset_index(drop=True)
        if df_gwas_sub.shape[0] == 0:
            continue
        
        logging.info(f'Chromosome {i}: Loading genotype covariance meta information.')
        df_cov_meta = load_cov_meta(args.genotype_covariance.format(chr_num=i))
        
        # step0
        n0 = df_weight_sub.shape[0]  # for book keeping
        # we enforce the GWAS table and the IDP weights to have 
        # the same SNPs as genotype covariance
        # the weights of the missing ones are set to NaN.
        df_gwas_sub = rearrage_df_by_target(
            df=df_gwas_sub, 
            target=df_cov_meta,
            df_value_cols=['effect_size']
        )
        df_weight_sub = rearrage_df_by_target(
            df=df_weight_sub, 
            target=df_cov_meta,
            df_value_cols=list(df_weight.columns[4:])
        )
        n1 = df_gwas_sub.effect_size.notna().sum()
        logging.info('Step0 Chromosome {}: {} out of {} SNPs in IDP/GWAS are used.'.format(i, n1, n0))
        
        logging.info(f'Step1 Chromosome {i}: Working with genotype covariance.')
        
        weight = df_weight_sub.iloc[:, 4 : ].to_numpy(copy=True)
        weight[np.isnan(weight)] = 0
        
        b_gwas = df_gwas_sub.effect_size.to_numpy(copy=True)
        b_gwas[np.isnan(b_gwas)] = 0
        se_gwas = df_gwas_sub.effect_size_se.to_numpy(copy=True)
        se_gwas[np.isnan(se_gwas)] = 1
        z_gwas = b_gwas / se_gwas
        
        cov_mat = CovMatrix(args.genotype_covariance.format(chr_num=i))
        cov_x_weight, diag_cov = cov_mat.eval_matmul_on_left(weight)
        D_chr = weight.T @ cov_x_weight
        del cov_x_weight
        var_R_chr = diag_cov
        numer_b_chr = weight.T @ ( var_R_chr * b_gwas )
        numer_z_chr = weight.T @ ( np.sqrt(var_R_chr) * z_gwas )
        D += D_chr
        numer_b += numer_b_chr
        numer_z += numer_z_chr
        
    logging.info('Step2: Computing marginal test.')
    S_D = np.sqrt(D.diagonal())
    beta_imagexcan = numer_b / np.power(S_D, 2)
    z_imagexcan = numer_z / S_D
   
    logging.info('Step3: Running susieR.')
    Sigma =  D / S_D[:, np.newaxis] / S_D[np.newaxis, ]
    susie_pip, susie_cs = run_susie_wrapper(z_imagexcan, Sigma, params={'z_ld_weight': args.z_ld_weight})
         
    logging.info('Saving outputs.')
    df_res = pd.DataFrame({
        'IDP': idp_names,
        'bhat': beta_imagexcan,
        'pval': z2p(z_imagexcan),
        'pip': susie_pip,
        'cs95': susie_cs
    })
    df_res.to_csv(args.output, index=False)
    
    logging.info('Done.')
    
