import pandas as pd
from run_gw_ridge import load_genotype_from_bedfile

def add_noise(y, sd_noise):
    return y + np.random.normal(scale=sd_noise, size=(y.shape[0]))

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='simulate_phenotypes.py', description='''
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
    parser.add_argument('--rand_seed', nargs='+', type=int, help='''
        The list of random seeds  
    ''')
    parser.add_argument('--num_mediators', default=None, type=int, help='''
        Number of mediators  
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
    parser.add_argument('--pi0', default=None, type=float, help='''
        pi0 in [0, 1]. The probability of beta_k being zero.
    ''')
    parser.add_argument('--param_config', default=None, help='''
        Parameters in one config file. Should include pi0, h2s, pves, num_mediators
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
        h2s, pves, pi0, num_mediators = args.h2s, args.pves, args.pi0, args.num_mediators
    else:
        with open(args.param_config, 'r') as f:
            try:
                param_dict = yaml.safe_load(f)
            except:
                raise ValueError('Something wrong in param_config')
        h2s = param_dict['h2s'] if args.h2s is None
        pves = param_dict['pves'] if args.pves is None
        pi0 = param_dict['pi0'] if args.pi0 is None
        num_mediators = param_dict['num_mediators'] if args.num_mediators is None
        
    logging.info('Loading meta information')
    # load individual list
    samples = load_indiv(args.indiv_list)
    # get genotype meta data
    snpids = get_geno_meta(args.geno_bed_pattern)
    nsnp = len(snpids)
    
    logging.info('Simulation effect sizes')
    B = np.random.normal(size=(nsnp, num_mediators))
    beta = np.random.normal(size=(num_mediators))
    beta[np.random.uniform() < pi0] = 0
    b_null = np.random.normal(size=nsnp)
    df_mediator = pd.DataFrame({'mediator': [ i for i in range(num_mediators) ], 'beta': beta})
    
    logging.info('Calculating genetic component of mediators/y_null')
    start = 0
    gm = np.zeros((len(indiv_list), num_mediators))
    ynullm = np.zeros(len(indiv_list))
    df_snp = []
    for i in range(1, 23):
        geno_prefix = geno_bed_pattern.format(chr_num=i)
        geno_i, _, _, snp_meta = load_genotype_from_bedfile(
            f'{geno_prefix}.bed', indiv_list, snplist_to_exclude=set([]), 
            return_snp=True, standardize=True)
        df_dic = {
            'snpid': snp_meta[0], 'chr': snp_meta[3], 
            'a0': snp_meta[1], 'a1': snp_meta[2]}
        B_i = B[start : (start + geno_i.shape[1]), :]
        b_i = b_null[start : (start + geno_i.shape[1])]
        gm += geno_i @ B_i
        ynullm += geno_i @ b_i
        for k in range(num_mediators):
            df_dic[f'B_{k}'] = B_i[:, k]
        df_dic['b_y_null'] = b_null
        df_snp.append(pd.DataFrame(df_dic))
        start += geno_i.shape[1]
    df_snp = pd.concat(df_snp, axis=0)
    
    df_gmed = pd.concat([ 
        pd.DataFrame({'individual': indiv_list}),
        pd.DataFrame(gm, columns=[ f'm_{k}' for k in range(num_mediators)],
        axis=1)
    df_gynull = pd.DataFrame({'individual': indiv_list, 'gynull': ynullm})
    
    logging.info('Simulating observed mediators/y_null')
    var_gmed = df_gmed.iloc[:, 1:].var(axis=0)
    gmed_names = var_gmed.index.tolist()
    gmed_values = var.gmed.values.tolist()
    df_omed = {'individual': indiv_list}
    for h2 in h2s:
        for n, v in zip(gmed_names, gmed_values):
            error_sd = np.sqrt(v / h2 * (1 - h2))
            df_omed[f'{n}_h2_{h2}'] = add_noise(df_gmed[[n]].values, error_sd)        
    df_omed = pd.DataFrame(df_omed)
    df_y = {'individual': indiv_list}
    gynull = df_gynull[['gynull']].values
    for h2 in h2s:
        for pve in pves:
            ff = h2 * pve
            error_sd = np.sqrt(gynull.var() / ff * (1 - ff))
            df_y[f'null_h2_{h2}_pve_{pve}'] = add_noise(gynull, error_sd)
    # df_oynull = pd.DataFrame(df_oynull)
    
    logging.info('Simulating y')
    # df_oy = {'individual': indiv_list}
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
    df_snp.to_parquet(f'{args.output_prefix}.rand_{args.rand_seed}.snp_effect.parquet')
    df_mediator.to_parquet(f'{args.output_prefix}.rand_{args.rand_seed}.mediator_effect.parquet')
    df_gmed.to_parquet(f'{args.output_prefix}.rand_{args.rand_seed}.gmed.parquet')
    df_omed.to_parquet(f'{args.output_prefix}.rand_{args.rand_seed}.omed.parquet')
    df_y.to_parquet(f'{args.output_prefix}.rand_{args.rand_seed}.oy.parquet')
    
    
        
            
    