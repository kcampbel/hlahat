#!/usr/bin/env python
""" Timestamp expression matrix, metadata and update log
""" 
import os
import sys
import pandas as pd 
import argparse
import logging
import time
import re
from process_exprs.fileio import read_exprs, write_exprs, get_extension, file_time, append_basename

def parse_args(args=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument('--exprs', required=True, help='Gene expression matrix file')
    parser.add_argument('--exprs_log', required=True, help='Log file')
    parser.add_argument('--metadata_tsv', help='Metadata file')
    parser.add_argument('-o', '--outpath', required=True, help='Output path (e.g. s3://)')

    return parser.parse_args(args)

def rev_filename(filename, timestamp, outpath):
    fn = append_basename(filename, timestamp)
    fp = f'{outpath}/timestamp/{fn}'
    fp_orig = f'{outpath}/{os.path.basename(filename)}'
    return(fp_orig, fp)

def main():
    formatter = '%(asctime)s:%(levelname)s:%(name)s:%(funcName)s: %(message)s'
    logging.basicConfig(format=formatter, level=logging.INFO)
    logging.info(f'Starting {os.path.basename(__file__)}')

    args = parse_args()
#    args = parse_args(
#        [
#        '--exprs', '/home/csmith/git/bioinfo-fio/tme/test/eset_pact_geneid_proteincoding.tsv',
#        '--exprs_log', '/home/csmith/git/bioinfo-fio/tme/test/pact_eset_log.tsv',
#        '--metadata_tsv', '/home/csmith/git/bioinfo-fio/tme/test/meta_eset.tsv',
#        #'--metadata_tsv', '/home/csmith/git/bioinfo-fio/tme/test/meta_eset_20210625T165418.tsv', # Has PACT004
#        #'--exprs', '/home/csmith/git/bioinfo-fio/tme/test/eset_pact_geneid_proteincoding_20210625T165418.tsv', # Has PACT004
#        '-o', '/tmp/update_exprs',
#        ]
#    )

    # S3
#    args = parse_args(
#        [
#        '--exprs', '/home/csmith/git/bioinfo-fio/tme/test/eset_pact_geneid_proteincoding.tsv',
#        '--exprs_log', '/home/csmith/git/bioinfo-fio/tme/test/pact_eset_log.tsv',
#        '--metadata_tsv', '/home/csmith/git/bioinfo-fio/tme/test/meta_eset.tsv',
#        #'--metadata_tsv', '/home/csmith/git/bioinfo-fio/tme/test/meta_eset_20210625T165418.tsv', # Has PACT004
#        #'--exprs', '/home/csmith/git/bioinfo-fio/tme/test/eset_pact_geneid_proteincoding_20210625T165418.tsv', # Has PACT004
#        '-o', 's3://pact-research/csmith/tme/test/timestamp_exprs',
#        ]
#    )
    # Expression set
    exprs = read_exprs(args.exprs, index_col=0)
    timestamp = file_time(args.exprs)
    fp_orig, fp = rev_filename(args.exprs, timestamp, args.outpath)
    logging.info(f'Writing {fp} and {fp_orig}')
    dn = os.path.dirname(fp)
    if not os.path.exists(dn) and not re.match('s3://', fp):
        os.makedirs(dn, exist_ok=True)
    write_exprs(exprs, fp, compression='gzip')
    write_exprs(exprs, fp_orig, compression='gzip')

    if args.metadata_tsv:
        meta = pd.read_csv(args.metadata_tsv, sep='\t', index_col=0)
        fp_orig, fp = rev_filename(args.metadata_tsv, timestamp, args.outpath)
        logging.info(f'Writing {fp} and {fp_orig}')
        meta.to_csv(fp, sep='\t')
        meta.to_csv(fp_orig, sep='\t')

    # Log
    exprs_log = pd.read_csv(args.exprs_log, sep='\t')
    fp = f'{args.outpath}/{os.path.basename(args.exprs_log)}'
    logging.info(f'Updating {fp}')
    exprs_log.to_csv(fp, sep='\t', index=False)

    logging.info(f'{os.path.basename(__file__)} finished')


if __name__ == '__main__':
    main()
