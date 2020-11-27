import pandas as pd
import numpy as np


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
    

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='run_simagexcan.py', description='''
        Run S-ImageXcan with pre-computed genotype covariance.
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
    parser.add_argument('--output_prefix', help='''
        The prefix of output.
        Will return both marginal test result and also the susieR result.
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
    
    logging.info('Loading GWAS.')
    df_gwas = load_gwas(args.gwas)
    # df_gwas columns: 
    # snpid, non_effect_allele, effect_allele, 
    # effect_size, effect_size_se, chr
    logging.info('GWAS SNP = {}'.format(df_gwas.shape[0]))
    
    logging.info('Loading IDP weights.')
    df_weight = load_idp(args.idp_weight)
    logging.info('IDP SNP = {}'.format(df_weight.shape[0]))
    
    logging.info('Harmonizing GWAS and IDP weights.')
    # harmonize GWAS and IDP weight table so that they have the same set of 
    # SNPs (including direction).
    df_gwas, df_weight = harmonize_gwas_and_weight(df_gwas, df_weight)
    logging.info('{} SNPs left after harmonizing GWAS and IDP weights.'.format(df_gwas.shape[0]))
    
    for i in tqdm(range(1, 23)):
    
        logging.info(f'Chromosome {i}: Loading genotype covariance meta information.')
        df_cov_meta = load_cov_meta(args.genotype_covariance.format(chr_num=i))
    
        # please refer to https://github.com/hakyimlab/yanyu-notebook/blob/master/notes/date_112420.Rmd
        # for the details of the S-ImageXcan formula
        # to take the following procedure.
        # 0. subset weight to GWAS SNPs.
        # 1. obtain D and S_R.
        # 2. compute marginal test.
        # 3. run susieR.
        
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
        logging.info('Chromosome {i} --> Step0: {} out of {} SNPs in IDP/GWAS are used.'.format(n1, n0)
        
        logging.info('Chromosome {i} --> Step1: Working with genotype covariance.')
        
    
    
    
    
    
