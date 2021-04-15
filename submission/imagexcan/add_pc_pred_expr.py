if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='add_pc_weights.py', description='''
        Add PC weights in input_pc to current input data.frame
    ''')
    parser.add_argument('--input', help='''
        PARQUET
    ''')
    parser.add_argument('--input_pc', help='''
        PARQUET
    ''')
    parser.add_argument('--output', help='''
        PARQUET
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

    logging.info('Loading input.')
    df1 = pd.read_parquet(args.input)
    logging.info('input shape = {}'.format(df1.shape))
    
    logging.info('Loading input_pc.')
    df2 = pd.read_parquet(args.input_pc)
    logging.info('input shape = {}'.format(df2.shape))
    
    logging.info('Getting PC columns.')
    pc_cols = []
    for i in df2.columns:
        if 'PC-' in i:
            pc_cols.append(i)
    logging.info('There are {} PCs.'.format(len(pc_cols)))
    df2 = df2[['indiv'] + pc_cols].copy()
    
    logging.info('Adding input_pc to input.')
    df_merge = pd.merge(df1, df2, on=['indiv'], how='inner')
    logging.info('There are {} NAs in the resulting data frame'.format(df_merge.isna().sum().sum())
    if df_merge.indiv.unique().shape[0] != df_merge.shape[0]:
        raise ValueError('Duplicated indiv.')
    logging.info('Merged df shape = {}'.format(df_merge.shape))
    
    logging.info('Saving to output.')
    df_merge.to_parquet(args.output)
    
    
    
