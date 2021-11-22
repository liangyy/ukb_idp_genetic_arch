import pandas as pd
import numpy as np

def load_b(args_):
    df_b = pd.read_parquet(args_[0])
    for i in args_[1:]:
        tmp = i.split(':')
        df_b.rename(columns={tmp[1]: tmp[0]}, inplace=True)
    return df_b
def load_list(fn):
    res = []
    with open(fn, 'r') as f:
        for i in f:
            i = i.strip()
            res.append(i)
    return res

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='calc_assoc.py', description='''
        Calculate b1' X'X b2 and b2' X'X b2 with standarized genotype X. 
    ''')
    parser.add_argument('--genotype_covariance', help='''
        The genotype covariance computed in build_genotype_covariance.py
        Accept wildcard {chr_num}.
        Will automatically search for the corresponding meta SNP file.
    ''')
    parser.add_argument('--b1', nargs='+', help='''
        Need to have column names for: 
            snpid, non_effect_allele, effect_allele, chr.
        The format is: snpid:rsid_col, ...
    ''')
    parser.add_argument('--b1_effs', help='''
        The list of effect size columns in --b1
    ''')
    parser.add_argument('--b2', nargs='+', help='''
        Need to have column names for: 
            snpid, non_effect_allele, effect_allele, chr.
        The format is: snpid:rsid_col, ...
    ''')
    parser.add_argument('--b2_effs', help='''
        The list of effect size columns in --b2
    ''')
    parser.add_argument('--no_standarize', action='store_true', help='''
        If specified, it will skip genotype standardization.
    ''')
    parser.add_argument('--output', help='''
        The output CSV filename.
        Will return b1' X'X b2 and b2' X'X b2 for all b1_effs and b2_effs
        combinations 
    ''')
    args = parser.parse_args()
    
    import logging, time, sys, os
    # configing util
    logging.basicConfig(
        level = logging.INFO, 
        stream = sys.stderr, 
        format = '%(asctime)s  %(message)s',
        datefmt = '%Y-%m-%d %I:%M:%S %p'
    )
    from CovConstructor import CovMatrix
    from run_simagexcan import load_cov_meta, rearrage_df_by_target
    
    b1_effs = load_list(args.b1_effs)
    b2_effs = load_list(args.b2_effs)
    
    logging.info('Loading b1 and b2')
    df_b1 = load_b(args.b1)
    df_b2 = load_b(args.b2)
    n1 = df_b1.shape[0]
    n2 = df_b2.shape[0]
    
    logging.info('Calculating by chromosome')
    n1_used = 0
    n2_used = 0
    res12 = np.zeros((len(b1_effs), len(b2_effs)))
    res22 = np.zeros((len(b2_effs)))
    for i in range(1, 23):
        logging.info(f'Working on chromosome {i}')
        df_meta = load_cov_meta(args.genotype_covariance.format(chr_num=i))
        df_b1_sub = pd.merge(df_meta[['snpid', 'chr']], df_b1, on='snpid')
        df_b2_sub = pd.merge(df_meta[['snpid', 'chr']], df_b2, on='snpid')
        breakpoint()
        df_b1_sub = rearrage_df_by_target(
            df=df_b1_sub,
            target=df_meta,
            df_value_cols=b1_effs)
        df_b2_sub = rearrage_df_by_target(
            df=df_b2_sub,
            target=df_meta,
            df_value_cols=b2_effs)
        n1_used += df_b1_sub[b1_effs[0]].notna().sum()
        n2_used += df_b2_sub[b2_effs[0]].notna().sum()
        
        cov_mat = CovMatrix(args.genotype_covariance.format(chr_num=i))
        
        b1_values = df_b1_sub[b1_effs].values.copy()
        b2_values = df_b2_sub[b2_effs].values.copy()
        if not args.no_standarize:
            logging.info('Standardizing genotypes')
            _, diag_cov = cov_mat.eval_matmul_on_left(b2_values, param=100)
            b2_values = b2_values / diag_cov.sqrt()[:, np.newaxis]
            b1_values = b1_values / diag_cov.sqrt()[:, np.newaxis]
        cov_x_b2, diag_cov = cov_mat.eval_matmul_on_left(b2_values, param=100)
        b1_cov_x_b2 = b1_values.T @ cov_x_b2
        b2_cov_x_b2 = np.einsum(b2_values, cov_x_b2, 'ij,ij->j')
        res12 += b1_cov_x_b2
        res22 += b2_cov_x_b2    
    logging.info(f'{n1_used} out of {n1} SNPs in b1 are used')
    logging.info(f'{n2_used} out of {n2} SNPs in b2 are used')
    
    logging.info('Saving as data frames')
    df = []
    for i, name in enumerate(b1_effs):
        tmp = pd.DataFrame({
            'left': [ name for j in range(res12.shape[1]) ],
            'right': b2_effs,
            'left_cov_x_right': res12[i, :]})
        df.append(tmp)
    df.append(
        pd.DataFrame({
            'left': b2_effs,
            'right': b2_effs,
            'left_cov_x_right': res22}))
    df = pd.concat(df, axis=0)
    df.to_csv(args.output, index=False)
    
