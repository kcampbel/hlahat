#!/usr/bin/env python
""" Retrieve TME reference data """
import os
import argparse
import re
import subprocess as sb
import yaml
import numpy as np
import logging
from gtfparse import read_gtf
from tme import data
from tme.fileio import package_file_path, get_file

def parse_args(args=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument('--install_config', default=package_file_path(data, 'install.yml'), help='References yml')
    #parser.add_argument('--params', default=package_file_path(data, 'tme_report_params.yml'), help='TME parameters yaml')
    parser.add_argument('--outpath', default=package_file_path(data, ""), help='GTF file')
    
    return parser.parse_args(args)

def main():
    logging.info(f'Starting {os.path.basename(__file__)}')
    args = parse_args()
#    args = parse_args(
#        [
#        '--config', 'test_data/get_reference.yml',
#        '--outpath', '/tmp'
#        ]
#    )
    logging.info(args)

    ref = dict()
    install_config = yaml.safe_load(open(args.install_config))
    # Disable writing to params file for now
    #params = yaml.safe_load(open(args.params))

    # GTF
    cols = ['gene_name', 'gene_id', 'seqname', 'start', 'end', 'feature', 'source', 'gene_biotype']
    gtf_f = get_file(install_config['gtf'], '/tmp')
    gtf = read_gtf(gtf_f, usecols=cols)
    gtf_out = gtf[
        (gtf.feature == 'gene') &
        (gtf.gene_biotype == 'protein_coding') &
        (gtf.source.str.contains('ensembl')) &
        (~gtf.seqname.str.contains('GL'))
    ]\
      .sort_values('gene_name')[cols]\
      .rename(columns={'seqname': 'chromosome'})

    gtf_out['chromosome'] = np.where(
        gtf_out.chromosome == 'MT', 
        'chrM', 'chr' + gtf_out.chromosome
    )
    fn = f'{args.outpath}/{os.path.basename(gtf_f).replace("gtf", "tsv")}'
    logging.info(f'Writing {fn}')
    gtf_out.to_csv(fn, sep='\t', index=False)
    #params['gtf_tsv_f'] = fn

    # HGNC
    hgnc_f = get_file(install_config['hgnc'], args.outpath)
    fn = f'{args.outpath}/{os.path.basename(hgnc_f)}'
    if not fn.endswith('.gz'):
        sb.run(['gzip', '-f', fn], check=True)
        fn = fn + '.gz'
    
    # Write yaml with reference locations
    #params['hgnc_f'] = fn
    ##fn = f'{package_file_path(data, "")}/reference.yml'
    #logging.info(f'Writing {args.params}')
    #with open(args.params, 'wt') as fh:
    #    yaml.dump(params, fh)

if __name__ == '__main__':
    main()
