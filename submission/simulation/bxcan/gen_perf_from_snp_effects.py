if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='gen_perf_from_snp_effects.py', description='''
        Generate a fake performance table for a SNP effect table so that 
        we can run BrainXcan using these SNP effects.
    ''')
    parser.add_argument('--parquet', help='''
        SNP effect parquet.
    ''')
    parser.add_argument('--output', help='''
        Output TSV.GZ performance file.
    ''')
    parser.add_argument('--output', help='''
        Output file (TAB-delimited) including Tr(R'R) and Tr^2(R) for each 
        chromosome. 
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
    
    logging.info('Loading parquet')
    df = pd.read_parquet(args.parquet)
    phenos = list(df.columns)[4:]
    
    logging.info('Generating fake performance table')
    df_perf = pd.DataFrame({
        'R2': [ 0 for i in phenos ],
        'Pearson': [ 0 for i in phenos ],
        'Spearman': [ 0 for i in phenos ],
        'phenotype': phenos})
    df_perf.to_csv(args.output, compression='gzip', sep='\t', index=False)
    logging.info('Done')
    