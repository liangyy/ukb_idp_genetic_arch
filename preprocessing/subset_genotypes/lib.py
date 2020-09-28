def get_bed_files(config):
    if 'genotype_bed' not in config:
        raise ValueError('No genotype_bed in config.')
    file_prefix = config['genotype_bed']
    files = [ file_prefix + '.' + ss for ss in ['bim', 'bed', 'fam'] ]
    command = '--bfile ' + file_prefix
    return files, command