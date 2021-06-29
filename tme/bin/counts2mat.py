#!/usr/bin/env python
""" Counts to expression matrix converter.
Convert long format tsvs (e.g. RSEM) to a matrix with gene ids in rows and sample names in columns.
File name must be {sample_name}_tumor to parse the sample name for the column.
e.g. for RNA_PACT999_T_9999_tumor_rna.tar.gz, PACT999_T_9999 will be used as the column name.
""" 
import os
import sys
import pandas as pd 
import argparse
import logging
from process import counts2mat

def parse_args(args=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument('metric', choices=['TPM', 'FPKM', 'counts'], help='Count metric')
    parser.add_argument('counts', nargs='+', help='count files to process (e.g. *genes.tsv')
    parser.add_argument('-o', '--outfile', required=True, help='Output file')
    parser.add_argument('-u', '--update', type=str, help='Existing expression matrix file name to update')
    parser.add_argument('-n', '--gene_name', default='hgnc_symbol', help='Gene name column in counts file')
    parser.add_argument('-f', '--force', action='store_true', help='Overwrite existing samples in eset')
    parser.add_argument('--fof', action='store_true', help='counts is a file of files')
    
    return parser.parse_args(args)

def main():
    args = parse_args()
#    args = parse_args(
#        [
#            'TPM', 
#            './test/RNA_PACT004_T_560351F_tumor_rna.genes.tsv', 
#            '--gene_name', 'gene_id',
#            #'--fof',
#            #'-o', './test/eset_pact_geneid_newsample.tsv',
#            '-u', './test/eset_pact_geneid_newsample.tsv',
#            #'-u', './test/eset_pact_geneid.tsv',
#            '-o', '/dev/null',
#            '-f',
#        ])
    formatter = '%(asctime)s:%(levelname)s:%(name)s:%(funcName)s: %(message)s'
    logging.basicConfig(format=formatter, level=logging.DEBUG)
    logging.info(f'Starting {os.path.basename(__file__)}')

    # Collect new counts to add to eset
    if args.fof:
        reader = open(args.counts[0]).readlines()
        counts = [x.strip('\n') for x in reader]
    else:
        counts = args.counts
    
    # Load eset if updating
    if args.update:
        eset = pd.read_csv(args.update, sep='\t', index_col=0)
    else:
        eset = pd.DataFrame()
        eset.index.name = args.gene_name

    df = counts2mat(counts, args.metric, args.gene_name, eset, args.force)
    if df is not None:
        logging.info(f'Writing {args.outfile}')
        if 'parquet' in args.outfile:
            df.to_parquet(args.outfile)
        else:
            df.to_csv(args.outfile, sep='\t')
    logging.info(f'{os.path.basename(__file__)} finished')

if __name__ == '__main__':
    main()
