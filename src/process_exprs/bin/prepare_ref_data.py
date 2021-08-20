#!/usr/bin/env python
""" Process_exprs reference data 

Reference data:
- gene_id -> symbol mapping file for Xena emat
"""
import os
import argparse
import pandas as pd
import yaml
import numpy as np
import logging
from gtfparse import read_gtf
from process_exprs import data
from commonLib.lib.fileio import package_file_path, get_file

def parse_args(args=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument('--config', default=package_file_path(data, 'pipeline.yml'), help='References yml')
    parser.add_argument('--outfile', default=f'{package_file_path(data, "")}/hgnc_pact_xena.tsv', help='Output filename')
    parser.add_argument('--force', action='store_true', help='Force overwrite')
    
    return parser.parse_args(args)

def munge_gtf(df, cols:list):
    """ Munge a GTF and return cols """
    out = df[
        (df.feature == 'gene') 
    ]\
      .sort_values('gene_name')\
      .rename(columns={'seqname': 'chromosome'})\

    out['chromosome'] = np.where(
        out.chromosome == 'MT', 
        'chrM', 'chr' + out.chromosome
    )
    out = out[cols].drop_duplicates()
    return out


def main():
    logging.info(f'Starting {os.path.basename(__file__)}')
    args = parse_args()
    logging.info(args)
    if os.path.exists(args.outfile) and not args.force:
        raise FileExistsError(f'{args.outfile} exists. Set --force to overwrite')

    config = yaml.safe_load(open(args.config))
    config = config['install']

    # Generate a gene id -> hgnc mapping file for Xena emat
    ## GTF
    gtf_f = f'/tmp/{config["gtf"]}'
    if os.path.exists(gtf_f):
        logging.info(f'{gtf_f} exists. Using existing file.')
    else:
        gtf_f = get_file(config['gtf'], '/tmp')

    gtf_df = read_gtf(gtf_f)
    hgnc_pact = munge_gtf(gtf_df, cols=['gene_id', 'gene_name'])

    ## Xena Gene ID -> Symbol mapping
    xena_map_f = 'https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/probeMap%2Fgencode.v23.annotation.gene.probemap'
    hgnc_xena = pd.read_table(xena_map_f)
    hgnc_xena['gene_id'] = hgnc_xena['id'].str.split('.').str[0]
    hgnc_xena = hgnc_xena.rename(columns={'gene':'gene_name'})\
        .loc[:, ['gene_id', 'gene_name']]\
        .drop_duplicates()
    
    ## Merge
    logging.info('Merging PACT and Xena gene_id -> gene_name mappings')
    hgnc_m = hgnc_pact.merge(hgnc_xena, how='outer', indicator=True)
    hgnc_m['_merge2'] = np.select(
        [
            hgnc_m._merge == 'left_only',
            hgnc_m._merge == 'right_only'
        ],
        ['pact_only', 'xena_only'], hgnc_m._merge
    )
    logging.info(f'Gene id merge: \n{hgnc_m._merge2.value_counts()}')

    logging.info('Merge by retaining mapping in both Xena and PACT + those unique to PACT. Then add mappings unique to Xena')
    lo = hgnc_m[(hgnc_m._merge == 'left_only') | (hgnc_m._merge == 'both')]
    ro = hgnc_m[hgnc_m._merge == 'right_only']
    ro = ro[~ro.gene_id.isin(lo.gene_id)]

    loro = pd.concat([lo, ro], ignore_index=True)
    logging.info(f'Gene id merge: \n{loro._merge2.value_counts()}')

    logging.info('Find remaining gene id dups')
    dups = loro.gene_id.value_counts().reset_index(name='cnt')
    dups = dups[dups.cnt>1]
    logging.info(f'Remaining dups that need to be manually curated: \n{dups}')

    logging.info(f'Writing {args.outfile}')
    loro.to_csv(args.outfile, sep='\t', index=False)
    logging.info(f'{os.path.basename(__file__)} finished')

if __name__ == '__main__':
    main()
