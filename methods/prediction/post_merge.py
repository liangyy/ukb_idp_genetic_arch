import pandas as pd

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='post_merge.py', description='''
        This is a standalone script merging the predicted idps from 1 .. 22 
        chromosomes into one parquet file.
    ''')
    parser.add_argument('--input_pattern', help='''
        Contain {chr_num} as wildcard.
        Will load chr_num = 1 to 22.
    ''')
    parser.add_argument('--indiv_col', default='indiv', help='''
        Column name of individual ID. 
    ''')
    parser.add_argument('--output', help='''
        The output is saved as one parquet file.
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
    
    logging.info('On chromosome 1.')
    df = pd.read_parquet(args.input_pattern.format(chr_num=1))
    cols = df.columns
    for i in range(2, 23):
        logging.info(f'On chromosome {i}.')
        tmp = pd.read_parquet(args.input_pattern.format(chr_num=i))
        tmp2 = pd.DataFrame({args.indiv_col : df[args.indiv_col]})
        tmp = pd.merge(tmp2, tmp, on=args.indiv_col)
        tmp = tmp[cols]
        df.iloc[:, 1:] = df.iloc[:, 1:] + tmp.iloc[:, 1:]
    
    logging.info('Saving.')
    df.to_parquet(args.output)
    
