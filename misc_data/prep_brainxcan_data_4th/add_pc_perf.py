if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='add_pc_perf.py', description='''
        Add PC weights in input_pc to current input data.frame
    ''')
    parser.add_argument('--input', help='''
        TSV.GZ
    ''')
    parser.add_argument('--input_pc', help='''
        TSV.GZ
    ''')
    parser.add_argument('--output', help='''
        TSV.GZ
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
    df1 = pd.read_csv(args.input, compression='gzip', sep='\t')
    logging.info('input shape = {}'.format(df1.shape))
    
    logging.info('Loading input_pc.')
    df2 = pd.read_csv(args.input_pc, compression='gzip', sep='\t')
    logging.info('input shape = {}'.format(df2.shape))
    
    logging.info('Getting PC rows.')
    pc_rows = []
    for i in df2.phenotype:
        if 'PC-' in i:
            pc_rows.append(i)
    logging.info('There are {} PCs.'.format(len(pc_rows)))
    df2 = df2[ df2.phenotype.isin(pc_rows) ].reset_index(drop=True)
    
    logging.info('Adding input_pc to input.')
    df_merge = pd.concat([df1, df2], axis=0)
    df_merge.fillna('NA', inplace=True)
    logging.info('Merged df shape = {}'.format(df_merge.shape))
    
    logging.info('Saving to output.')
    df_merge.to_csv(args.output, compression='gzip', sep='\t', index=False)
    
    
    
