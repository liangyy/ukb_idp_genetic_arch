import yaml

def load_yaml(fn):
    with open(fn, 'r') as f:
        res = yaml.safe_load(stream)
    return res

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='generate_permuted_gwas.py', description='''
        Generate permuted GWAS by LD block
    ''')
    parser.add_argument('--config_yaml', help='''
        Config file for genetic correlation run
    ''')
    parser.add_argument('--gwas_name', help='''
        GWAS name being used
    ''')
    parser.add_argument('--nrepeat', type=int, default=10, help='''
        Number of repeats
    ''')
    parser.add_argument('--seed', type=int, default=1, help='''
        Random seed 
    ''')
    parser.add_argument('--chr_col', help='''
        Column name of chromosome
    ''')
    parser.add_argument('--pos_col', help='''
        Column name of position
    ''')
    parser.add_argument('--value_cols', nargs='+', help='''
        Column names of value columns (these will be permuted)
    ''')
    parser.add_argument('--new_config', help='''
        Path to the new config file
    ''')
    parser.add_argument('--gwas_outdir', 
        help='''
        Output directory of permuted GWAS
    ''')
    parser.add_argument('--ldblock', help='''
        LD block BED file (TAB-delimited and base0)
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
    import pandas as pd
    import numpy as np
    from brainxcan.bxcan.util.misc import file_exists
    from brainxcan.bxcan.run_sbrainxcan import load_ldblock, get_idxs_by_block
    
    logging.info('Loading config')
    config = load_yaml(args.config_yaml)
    gwas_fn = config['gwas_file_pattern'].format(gwas_name=args.gwas_name)
    gwas_basename = os.path.basename(config['gwas_file_pattern'])
    if not file_exists(gwas_fn):
        raise ValueError(f'GWAS file {gwas_fn} does not exist. Exit!')
    
    logging.info(f'Loading GWAS {gwas_fn}')
    df_gwas = pd.read_csv(gwas_fn, sep='\t', compression='gzip')
    df_meta = pd.DataFrame({
        'chr': [ re.sub('^chr', '', i) for i in df_gwas[args.chr_col] ],
        'pos': df_gwas[args.pos_col].astype(int),
        'raw_idx': [ i for i in range(df_gwas.shape[0]) ]})
    
    logging.info('Loading LD block file')
    df_ldblock = load_ldblock(args.ldblock)
    
    logging.info('Start permutation')
    np.random.seed(args.seed)
    raw_idxs_by_block = []
    df_values = df_gwas[args.value_cols]
    df_non_values = df_gwas.drop(columns=args.value_cols)
    for i in range(1, 23):
        df_meta_i = df_meta[ df_meta.chr == str(i) ].reset_index(drop=True)
        df_ldblock_i = df_ldblock[ df_ldblock.chr == str(i) ].reset_index(drop=True)
        snp_idxs = get_idxs_by_block(df_meta_i, df_ldblock_i)
        raw_idxs_by_block += [ list(df_meta.raw_idx.values[i]) for i in snp_idxs ]
    for i in range(1, args.nrepeat + 1):
        logging.info(f'-> Repeat = {i}')
        raw_idxs_permuted = []
        block_idxs = np.random.permutation(len(raw_idxs_by_block))
        for j in block_idxs:
            raw_idxs_permuted += raw_idxs_by_block[j]
        sorted_raw_idxs_permuted = list(np.unique(np.array(raw_idxs_permuted)))
        if len(sorted_raw_idxs_permuted) != len(raw_idxs_permuted):
            raise ValueError('Permuted raw idxs have duplicated values')
        df_new = pd.concat([
            df_non_values.iloc[sorted_raw_idxs_permuted, :], 
            df_values.iloc[raw_idxs_permuted, :]],
            axis=1)
        fn = gwas_basename.format(gwas_name=f'{args.gwas_name}_x_n{i}')
        logging.info(f'-> Saving to {fn}')
        df_new.to_csv(f'{args.gwas_outdir}/{fn}', index=False, sep='\t', compression='gzip')
    
    config['gwas_file_pattern'] = f'{args.gwas_outdir}/{gwas_basename}'
    logging.info('Generating new config YAML')
    with open(args.new_config, 'w') as f:
        yaml.dump(config, f, default_flow_style=False)
    
    logging.info('Done')
    