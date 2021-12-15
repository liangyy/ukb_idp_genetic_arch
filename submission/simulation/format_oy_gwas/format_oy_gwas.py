if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='format_oy_gwas.py', 
        description='''
            Format GWAS output from `gw_qtl` to txt.gz files
        ''')
    parser.add_argument('--pheno_list', help='''
        The list of phenotype file tag name to extract
    ''')
    parser.add_argument('--input_pattern', help='''
        The parquet file pattern with {chr_num} and {pheno} as wildcards
    ''')
    parser.add_argument('--snp_bim_pattern', help='''
        Corresponding SNP BIM file pattern with {chr_num} as wildcard
    ''')
    parser.add_argument('--output_pattern', help='''
        Output TXT.GZ file pattern with {pheno} as wildcard 
    ''')
    parser.add_argument('--sample_size', type=int, help='''
        Sample size to add as a column
    ''')
    args = parser.parse_args()
 
    import logging, time, sys, os
    # configing util
    logging.basicConfig(
        level = logging.INFO, 
        stream = sys.stderr, 
        format = '%(asctime)s  %(message)s',
        datefmt = '%Y-%m-%d %I:%M:%S %p')
    import pandas as pd
    
    logging.info('Loading phenotype list')
    pheno_list = []
    with open(args.pheno_list, 'r') as f:
        for i in f:
            i = i.strip()
            pheno_list.append(i)
    logging.info(f'There are {len(pheno_list)} phenotypes to work on')
    
    logging.info('Loading SNP meta info')
    df_snp = []
    for i in range(1, 23):
        fn = args.snp_bim_pattern.format(chr_num=i)
        tmp = pd.read_csv(fn, sep='\t', header=None)
        tmp.columns = ['chr', 'variant_id', 'ph', 'pos', 'alt', 'ref']
        df_snp.append(tmp)
    df_snp = pd.concat(df_snp, axis=0)
    
    logging.info('Loading GWAS')
    for pheno in pheno_list:
        logging.info(f'-> On phenotype = {pheno}')
        df_gwas = []
        for i in range(1, 23):
            fn = args.input_pattern.format(chr_num=i, pheno=pheno)
            tmp = pd.read_parquet(fn)
            tmp = pd.merge(
                tmp, 
                df_snp[['variant_id', 'chr', 'pos', 'alt', 'ref']], 
                on='variant_id')
            df_gwas.append(tmp)
        df_gwas = pd.concat(df_gwas, axis=0)
        df_gwas['sample_size'] = args.sample_size
        df_gwas.fillna('NA', inplace=True)
        df_gwas.to_csv(
            args.output_pattern.format(pheno=pheno), 
            compression='gzip',
            index=False,
            sep='\t')
    
    
