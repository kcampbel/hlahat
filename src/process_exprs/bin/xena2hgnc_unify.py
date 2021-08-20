#!/usr/bin/env python
""" Generate Xena HGNC emat 
Utilizes a PACT/Xena gene id -> symbol mapping file to generate an hgnc Xena expression matrix
"""
import pandas as pd
import logging
import os
import argparse
import re
from process_exprs import data
from commonLib.lib.fileio import package_file_path, get_file, read_exprs, write_exprs

def parse_args(args=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument('-i', '--infile', default='https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_rsem_gene_tpm.gz',
        help='Xena gene TPM input file')
    parser.add_argument('-o', '--outpath', required=True, help='Output path')
    parser.add_argument('-m', '--map', default=f'{package_file_path(data, "hgnc_pact_xena.tsv")}', help='PACT/Xena unified gene id to symbol mapping file')
    parser.add_argument('--force', action='store_true', help='Force overwrite')

    return parser.parse_args(args)    

def main():
    logging.info(f'Starting {os.path.basename(__file__)}')
    args = parse_args()
#    args = parse_args(
#        [
#           '-i', '/media/nfs/data/References/xena/TcgaTargetGtex_rsem_gene_tpm_one.gz',
#            '-o', '/tmp',
#        ]
#    )
    outfile = f'{args.outpath}/{re.sub(".gz", "_hgnc.parquet.gz", os.path.basename(args.infile))}'
    if os.path.exists(outfile) and not args.force:
        raise FileExistsError(f'{outfile} exists. Set --force to overwrite')

    hgnc = pd.read_csv(args.map,sep='\t', index_col='gene_id').loc[:, 'gene_name']

    logging.info(f'Reading {args.infile} and summarizing to gene from {args.map}')
    df = read_exprs(args.infile, index_col=0)
    df.index = df.index.str.rsplit('.').str[0]
    df = df.rpow(2) - 0.001
    dfm = df.merge(hgnc, how='left', left_index=True, right_index=True)\
        .groupby('gene_name').apply(sum, axis=0)
    dfm = dfm.drop(columns='gene_name')

    logging.info(f'Writing {outfile}')
    write_exprs(dfm, outfile, compression='gzip')
    logging.info(f'{os.path.basename(__file__)} finished')

if __name__ == '__main__':
    main()
