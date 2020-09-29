def get_bed_files(config):
    if 'genotype_bed' not in config:
        raise ValueError('No genotype_bed in config.')
    file_prefix = config['genotype_bed']
    files = [ file_prefix + '.' + ss for ss in ['bim', 'bed', 'fam'] ]
    command = '--bfile ' + file_prefix
    return files, command
    
def get_parquet_files(config):
    if 'genotype_parquet' not in config:
        raise ValueError('No genotype_parquet in config.')
    file_prefix = config['genotype_parquet']
    files = [ file_prefix + '.' + ss for ss in ['variants_metadata.parquet', 'variants.parquet'] ]
    command = file_prefix
    return files, prefix

def get_plink_filters(config):
    if 'plink_filters' not in config:
        return ''
    filters = config['plink_filters']
    return ' '.join(filters)
         