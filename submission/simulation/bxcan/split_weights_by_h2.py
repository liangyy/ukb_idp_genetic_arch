if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='split_weights_by_h2.py', description='''
        Split weights (output of train_ridge) by h2
    ''')
    parser.add_argument('--weight_prefix', help='''
        Prefix of weights (should have parquet and perf.tsv.gz files)
    ''')
    parser.add_argument('--h2s', type=str, nargs='+', help='''
        h2 being used
    ''')
    parser.add_argument('--output_prefix', help='''
        Output prefix
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

    logging.info(f'Loading weights = {args.weight_prefix}')
    weights = pd.read_parquet(args.weight_prefix + '.parquet')
    perf = pd.read_csv(
        args.weight_prefix + '.perf.tsv.gz', 
        sep='\t', compression='gzip')
    wnames = list(weights.columns.values[4:])
    wmeta = list(weights.columns.values[:4])
        
    for h2 in args.h2s:
        logging.info(f'Working on h2 = {h2}')
        h2_str = f'h2_{h2}'
        wnames_h2 = []
        for w in wnames:
            if h2_str in w:
                wnames_h2.append(w)
        logging.info(f'-> # weights = {len(wnames_h2)}')
        w_sub = weights[wmeta + wnames_h2]
        p_sub = perf[perf.phenotype.isin(wnames_h2)].reset_index(drop=True)
        logging.info(f'-> Saving w ({w_sub.shape[0]} x {w_sub.shape[1]}) and p ({p_sub.shape[0]} x {p_sub.shape[1]})')
        w_sub.to_parquet(
            f'{args.output_prefix}.{h2_str}.parquet', index=False)
        p_sub.to_csv(
            f'{args.output_prefix}.{h2_str}.perf.tsv.gz', compression='gzip', 
            sep='\t', index=False)
        

    
    
