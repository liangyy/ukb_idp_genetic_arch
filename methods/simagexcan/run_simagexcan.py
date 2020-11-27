import pandas as pd
import numpy as np
import pathlib
from collections import OrderedDict 
from pyutil import read_table


def rearrage_gwas_by_weight(df, target, df_value_cols):
    df_res = target[['snpid', 'chr', 'effect_allele', 'non_effect_allele']]
    df_res = pd.merge(
        df_res, df, 
        on=['snpid', 'chr'], 
        suffixes=['_res', '_df'],
        how='left'
    )
    flip_factor = check_flip(
        a1=df_target.effect_allele_res, 
        a2=df_target.non_effect_allele_res,
        b1=df_target.effect_allele_df, 
        b2=df_target.non_effect_allele_df
    )
    df_target[df_value_cols] = df_target[df_value_cols] * flip_factor[:, np.newaxis]
    df_target.drop(
        columns=['effect_allele_df', 'non_effect_allele_df'], inplace=True
    )
    df_target.rename(
        columns={
            'effect_allele_res': 'effect_allele',
            'non_effect_allele_res': 'non_effect_allele'
        }, 
        inplace=True
    )
    return df_target

def harmonize_gwas_and_weight(gwas, weight):
    '''
    Harmonize GWAS to weight SNP set.
    But only keep the ones that present in both.
    '''
    df_common = pd.merge(
        gwas[['snpid', 'chr', 'effect_allele', 'non_effect_allele']],
        weight[['snpid', 'chr', 'effect_allele', 'non_effect_allele']],
        on=['snpid', 'chr']
        suffixes=['_gwas', '_weight']
    )
    flip_factor = check_flip(
        a1=df_common.effect_allele_gwas, 
        a2=df_common.non_effect_allele_gwas,
        b1=df_common.effect_allele_weight, 
        b2=df_common.non_effect_allele_weight
    )
    df_gwas = pd.merge(
        df_common[['snpid', 'chr']], df_gwas, 
        on=['snpid', 'chr']
    )
    df_gwas.effect_size = df_gwas.effect_size * flip_factor
    df_weight = pd.merge(
        df_common[['snpid', 'chr']], df_weight, 
        on=['snpid', 'chr']
    )
    return df_gwas, df_weight

def _parse_gwas_args(args_list):
    fn = args_list[0]
    if not pathlib.Path(fn).is_file():
    raise ValueError('Filename is wrong. Cannot find the file.')
    dict = {}
    desired_cols = [
        'snpid', 'non_effect_allele', 'effect_allele', 
        'effect_size', 'effect_size_se', 'chr'
    ]
    snpid_name = None
    for i in args_list[1:]:
        
        tmp = i.split(':')
        if len(tmp) != 2:
            raise ValueError('Wrong gwas args list. Need [col]:[name] pairs.')
        if col not in desired_cols:
            raise ValueError(f'Wrong col = {col}.')
        col, name = tmp
        dict[col] = name
    rename_dict = OrderedDict()
    for dd in desired_cols:
        if dd not in dict:
            raise ValueError(f'Need to have col = {dd}.')
        rename_dict[dict[dd]] = dd
        if dd == 'snpid':
            snpid_name = dict[dd]
    return fn, rename_dict, snpid_name    
    
def load_gwas(gwas_args_list):
    fn, rename_dict, snpid_col = _parse_gwas_args(gwas_args_list)
    df = read_table(fn, indiv_col=snpid_col)
    df.rename(columns={'indiv_col': snpid_col}, inplace=True)
    df.rename(columns=rename_dict, inplace=True)
    return df[rename_dict.keys()]

def df_weight(fn):
    df = pd.read_parquet(fn)
    if 'a0' in df.columns:
        df.rename(columns={'a0': 'effect_allele'})
    if 'a1' in df.columns:
        df.rename(columns={'a1': 'non_effect_allele'})
    return df

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
    parser.add_argument('--gwas ', nargs='+', help='''
        Need to have column names for: 
            rsID, non_effect_allele, effect_allele, 
            effect_size, effect_size_se, chromosome.
        like: rsID:rsid_col, ..., chromosome:chr
    ''')
    parser.add_argument('--idp_weight', help='''
        The IDP weight table is in parquet format.
        It contains columns:
            snpid, effect_allele, non_effect_allele, chr.
        Along with all other columns for the IDPs.
    ''')
    parser.add_argument('--output', help='''
        The output CSV filename.
        Will return both marginal test result and also the susieR result.
    ''')
    parser.add_argument('--z_ld_weight', type=float, help='''
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
    from solver import run_susie_wrapper
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
    for i in tqdm(range(1, 23)):
    
        logging.info(f'Chromosome {i}: Loading genotype covariance meta information.')
        df_cov_meta = load_cov_meta(args.genotype_covariance.format(chr_num=i))
        
        # step0
        n0 = df_weight.shape[0]  # for book keeping
        # we enforce the GWAS table to have the same SNPs as the IDP weights
        # the weights of the missing ones are set to NaN.
        df_gwas = rearrage_df_by_target(
            df=df_gwas, 
            target=df_cov_meta
            df_value_cols=['effect_size']
        )
        df_weight = rearrage_df_by_target(
            df=df_gwas, 
            target=df_cov_meta
            df_value_cols=list(df_weight.columns[4:])
        )
        n1 = df_gwas.effect_size.notna().sum()
        logging.info('Step0 Chromosome {i}: {} out of {} SNPs in IDP/GWAS are used.'.format(n1, n0)
        
        logging.info('Step1 Chromosome {i}: Working with genotype covariance.')
        
        weight = df_weight.iloc[:, 4 : ].to_numpy(copy=True)
        weight[np.isnan(weight)] = 0
        
        b_gwas = df_gwas.effect_size.to_numpy(copy=True)
        b_gwas[np.isnan(b_gwas)] = 0
        se_gwas = df_gwas.effect_size.to_numpy(copy=True)
        se_gwas[np.isnan(se_gwas)] = 1
        z_gwas = b_gwas / se_gwas
        
        cov_mat = CovMatrix(args.genotype_covariance)
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
    beta_imagexcan = numer_b / np.pow(S_D, 2)
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
    
