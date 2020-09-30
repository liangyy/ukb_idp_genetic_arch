import pandas as pd

def read_list(filename):
    o = []
    with open(filename, 'r') as f:
        for i in f:
            o.append(i.strip())
    return o

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(prog='extract_from_genotype_parquet.py', description='''
        Extract a list of individuals and a list of SNPs 
        from genotype parquet file.
    ''')
    parser.add_argument('--input_prefix', help='''
        Prefix of input. 
        Will append '.variants.parquet' and '.variants_metadata.parquet'.
    ''')
    parser.add_argument('--snplist', help='''
        The list of SNP rsID to extract.
    ''')
    parser.add_argument('--indivlist', help='''
        The list of individual ID to extract.
    ''')
    parser.add_argument('--output_prefix', help='''
        Prefix of output genotype parquet.
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
   
    import pdb 

    logging.info('Loading genotype.')
    df_geno = pd.read_parquet(args.input_prefix + '.variants.parquet')
    
    logging.info('Loading genotype metadata.')
    df_meta = pd.read_parquet(args.input_prefix + '.variants_metadata.parquet')
    
    logging.info('Reading SNP list.')
    snplist = read_list(args.snplist)
    
    logging.info('Reading individual list.')
    indivlist = read_list(args.indivlist)
    
    logging.info('Subsetting individuals.')
    df_geno.drop(df_geno[ ~df_geno.individual.isin(indivlist) ].index, inplace=True)
    logging.info('Target individual list size = {}. Resulting genotype table has {} individuals remained.'.format(
        len(indivlist), df_geno.shape[0]
    ))
    
    logging.info('Subsetting SNPs in metadata table.')
    df_meta.drop(df_meta[ ~df_meta.rsid.isin(snplist) ].index, inplace=True)
    logging.info('Target SNP list size = {}. Resulting genotype table has {} SNPs remained.'.format(
        len(snplist), df_meta.shape[0]
    ))
    
    logging.info('Subsetting SNPs in genotype table.')
    df_geno.drop(columns=df_geno.columns.difference(['individual'] + df_meta.id.tolist()),inplace=True)
    
    logging.info('Writing output metadata table.')
    df_meta.to_parquet(args.output_prefix + '.variants_metadata.parquet')
    
    logging.info('Writing output genotype table.')
    df_geno.to_parquet(args.output_prefix + '.variants.parquet')
    
    logging.info('Done.')
    
