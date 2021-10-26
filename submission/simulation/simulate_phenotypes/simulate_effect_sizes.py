import pandas as pd
def get_geno_meta(pattern):
    n = 0
    for i in range(1, 23):
        bim = pattern.format(chr_num=i) + '.bim'
        snpids = []
        refs = []
        alts = []
        with open(bim, 'r') as f:
            for j in f:
                line = j.strip().split('\t')
                n, r, a = line[1], line[4], line[5]
                snpids.append(n)
                refs.append(r)
                alts.append(a)
    return pd.DataFrame({'snpid': snpids, 'ref': refs, 'alt': alts})
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='simulate_effect_sizes.py', 
        description='''
            Simulate effect size matrix B (normal), sparse effect size beta 
            (zero or normal). Also, we simulate the null by setting effect size 
            of SNPs on y as b (normal).
            Specify number of mediators to simulate.
            Set random seed.
            For each run, only one B, beta, and b are simulated.
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
    parser.add_argument('--rand_seed', type=int, help='''
        The list of random seeds  
    ''')
    parser.add_argument('--num_mediators', default=None, type=int, help='''
        Number of mediators  
    ''')
    parser.add_argument('--pi0', default=None, type=float, help='''
        pi0 in [0, 1]. The probability of beta_k being zero.
    ''')
    parser.add_argument('--param_config', default=None, help='''
        Parameters in one config file. Should include pi0, num_mediators
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
    if args.param_config is None:
        pi0 = args.pi0
        num_mediators = args.num_mediators
    else:
        with open(args.param_config, 'r') as f:
            try:
                param_dict = yaml.safe_load(f)
            except:
                raise ValueError('Something wrong in param_config')
        num_mediators = param_dict['num_mediators'] if args.num_mediators is None else args.num_mediators
        pi0 = param_dict['pi0'] if args.pi0 is None else args.pi0
    
    logging.info('Loading meta information')
    # load individual list
    indiv_list = load_indiv(args.indiv_list)
    # get genotype meta data
    df_snp = get_geno_meta(args.geno_bed_pattern)
 
    logging.info('Simulation effect sizes')
    B = np.random.normal(size=(nsnp, num_mediators))
    beta = np.random.normal(size=(num_mediators))
    beta[np.random.uniform(size=num_mediators) < pi0] = 0
    b_null = np.random.normal(size=nsnp)
    df_mediator = pd.DataFrame({'mediator': [ i for i in range(num_mediators) ], 'beta': beta})
    for k in range(num_mediators):
        df_snp[f'B_{k}'] = B[:, k]
    df_snp['b_y_null'] = b_null
    
    logging.info('Saving')
    df_snp.to_parquet(
        f'{args.output_prefix}.rand_{args.rand_seed}.snp_effect.parquet',
        index=False)
    df_mediator.to_parquet(
        f'{args.output_prefix}.rand_{args.rand_seed}.mediator_effect.parquet',
        index=False)
    