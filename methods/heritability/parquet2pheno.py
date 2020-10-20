def write_pheno(df, fn):
    with open(fn, 'w') as f:
        for i, p in zip(df.indiv.to_list(), df.pheno.to_list()):
            f.write(f'{i}\t{i}\t{p}\n')

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='parquet2pheno.py', description='''
        Read parquet file and generate the gcta pheno file.
    ''')
    parser.add_argument('--input', help='''
        Input parquet format.
    ''')
    parser.add_argument('--pheno_col', help='''
        Column name of the phenotype of interest.
    ''')
    parser.add_argument('--indiv_col', help='''
        Column name of the individual ID.
    ''')
    parser.add_argument('--output', help='''
        Output file name.
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
    
    logging.info('Loading phenotypes.')
    df = pd.read_parquet(args.input, columns=[args.indiv_col, args.pheno_col])
    
    df.rename(columns={args.indiv_col: 'indiv', args.pheno_col: 'pheno'}, inplace=True)
    
    logging.info('Writing to disk.')
    write_pheno(df, args.output)
   