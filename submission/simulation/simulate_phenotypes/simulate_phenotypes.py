import re
import numpy as np
import pandas as pd
from run_gw_ridge import load_genotype_from_bedfile

def add_noise(y, sd_noise):
    return y + np.random.normal(scale=sd_noise, size=(y.shape[0]))
def load_indiv(fn):
    res = []
    with open(fn, 'r') as f:
        for i in f:
            line = i.strip().split(' ')
            res.append(line[0])
    return res

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='simulate_phenotypes.py', description='''
        Given effect sizes,
        Simulate mediators with dense effect size matrix B (normal) and 
        phenotypes with sparse effect size beta (zero or normal).
        Also, we simulate the null by setting effect size of SNPs on y as 
        b (normal) with heritability = h2 * PVE.
        Specify h2, PVE, and number of mediators to simulate.
        Set random seed.
        For each run, only one B, beta, and b are simulated.
        Will simulate for each h2, PVE combination.
    ''')
    parser.add_argument('--geno_bed_pattern', help='''
        Genotype file in BED format (plink).
        It takes {chr_num} as wildcard.
        The script will load one chromosome at a time assuming the genotype file
        has 1 .. 22 chromosomes. 
    ''')
    parser.add_argument('--output_prefix', help='''
        Phenotype in parquet format.
    ''')
    parser.add_argument('--effect_size_prefix', help='''
        Prefix of the effect size parquet.
    ''')
    parser.add_argument('--rand_seed', type=int, help='''
        The list of random seeds  
    ''')
    parser.add_argument('--h2s', default=None, nargs='+', type=float, help='''
        The list of h2 
    ''')
    parser.add_argument('--pves', default=None, nargs='+', type=float, help='''
        The list of PVE
    ''')    
    parser.add_argument('--indiv_list', help='''
        The list of individuals
    ''')
    parser.add_argument('--param_config', default=None, help='''
        Parameters in one config file. Should include h2s, pves
    ''')
    args = parser.parse_args()
 
    import logging, time, sys, os
    # configing util
    logging.basicConfig(
        level = logging.INFO, 
        stream = sys.stderr, 
        format = '%(asctime)s  %(message)s',
        datefmt = '%Y-%m-%d %I:%M:%S %p')
    import yaml
    
    np.random.seed(args.rand_seed)
    
    # load parameters
    if args.param_config is None:
        h2s, pves = args.h2s, args.pves
    else:
        with open(args.param_config, 'r') as f:
            try:
                param_dict = yaml.safe_load(f)
            except:
                raise ValueError('Something wrong in param_config')
        h2s = param_dict['h2s'] if args.h2s is None else args.h2s
        pves = param_dict['pves'] if args.pves is None else args.pves
        
    logging.info('Loading effect sizes')
    # load individual list
    indiv_list = load_indiv(args.indiv_list)
    # load effect sizes
    df_snp = pd.read_parquet(args.effect_size_prefix + '.snp_effect.parquet')
    df_mediator = pd.read_parquet(args.effect_size_prefix + '.mediator_effect.parquet')
    num_mediators = df_mediator.shape[0]
    
    logging.info('Calculating genetic component of mediators/y_null')
    gm = np.zeros((len(indiv_list), num_mediators))
    ynullm = np.zeros(len(indiv_list))
    for i in range(1, 23):
        logging.info(f'-> Working on chromosome{i}')
        geno_prefix = args.geno_bed_pattern.format(chr_num=i)
        geno_i, _, _, snp_meta = load_genotype_from_bedfile(
            f'{geno_prefix}.bed', indiv_list, snplist_to_exclude=set([]), 
            return_snp=True, standardize=True)
        df_snp_i = pd.DataFrame({'snpid': snp_meta[0], 'a0': snp_meta[1], 'a1': snp_meta[2]})
        df_snp_i = pd.merge(
            df_snp_i, df_snp, 
            left_on=['snpid', 'a0', 'a1'], 
            right_on=['snpid', 'ref', 'alt'],
            how='left')
        df_snp_i.fillna(0, inplace=True)
        B_i = df_snp_i.iloc[:, 3:-1].values
        b_i = df_snp_i.iloc[:, -1].values
        gm += geno_i @ B_i
        ynullm += geno_i @ b_i
    
    df_gmed = pd.concat([ 
        pd.DataFrame({'individual': indiv_list}),
        pd.DataFrame(gm, columns=[ f'm_{k}' for k in range(num_mediators) ])],
        axis=1)
    df_gynull = pd.DataFrame({'individual': indiv_list, 'gynull': ynullm})
    
    logging.info('Simulating observed mediators/y_null')
    var_gmed = df_gmed.iloc[:, 1:].var(axis=0)
    gmed_names = var_gmed.index.tolist()
    gmed_values = var_gmed.values.tolist()
    df_omed = {'individual': indiv_list}
    for h2 in h2s:
        for n, v in zip(gmed_names, gmed_values):
            error_sd = np.sqrt(v / h2 * (1 - h2))
            df_omed[f'{n}_h2_{h2}'] = add_noise(df_gmed[n].values, error_sd)        
    df_omed = pd.DataFrame(df_omed)
    df_y = {'individual': indiv_list}
    gynull = df_gynull['gynull'].values
    for h2 in h2s:
        for pve in pves:
            ff = h2 * pve
            error_sd = np.sqrt(gynull.var() / ff * (1 - ff))
            df_y[f'null_h2_{h2}_pve_{pve}'] = add_noise(gynull, error_sd)
    
    logging.info('Simulating y')
    beta = df_mediator.iloc[:, 1].values
    for h2 in h2s:
        omed_cols = [ f'm_{k}_h2_{h2}' for k in range(num_mediators) ]
        omed_mat = df_omed[omed_cols].values
        my = omed_mat @ beta
        var_my = my.var()
        for pve in pves:
            error_sd = np.sqrt(var_my / pve * (1 - pve))
            oy = add_noise(my, error_sd)
            df_y[f'alt_h2_{h2}_pve_{pve}'] = oy
    df_y = pd.DataFrame(df_y)
    
    logging.info('Saving')
    df_gmed.to_parquet(
        f'{args.output_prefix}.rand_{args.rand_seed}.gmed.parquet',
        index=False)
    df_omed.to_parquet(
        f'{args.output_prefix}.rand_{args.rand_seed}.omed.parquet',
        index=False)
    df_y.to_parquet(
        f'{args.output_prefix}.rand_{args.rand_seed}.oy.parquet',
        index=False)
