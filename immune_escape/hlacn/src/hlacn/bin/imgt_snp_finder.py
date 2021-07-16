#!/usr/bin/env python
''' Generates SNPs from IMGT alignments
'''
import sys
import pandas as pd
import pickle
import argparse
import os
import logging
from hlacn.imgt import arr2snp, aln2snp, hla_nearest

def parse_args(args=None):
    parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='IMGT allele to SNP finder')
    parser.add_argument('alleles', nargs='+', help='Query alleles or file')
    parser.add_argument('-i', '--imgt_alignments', required=True,
        help='IMGT alignments pickle input')
    parser.add_argument('-g', '--ref_hla_genotypes', required=True, help='Reference HLA genotypes')
    parser.add_argument('-b', '--ref_hla_bed', required=True, help='Reference HLA bed')
    parser.add_argument('-o', '--outfile', default=sys.stdout)
    
    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)
#    args = parse_args(['--imgt_alignments', 'data/imgt_aln_gen.p',
#                       '-g', 'data/hla_genotypes_hs37d5.tsv',
#                       '-b', 'data/hla_hs37d5.bed',
#                       '-o', 'test/snps.tsv',
#                       'A*02:01']
#    )
    formatter = '%(asctime)s:%(levelname)s:%(name)s:%(funcName)s: %(message)s'
    logging.basicConfig(format=formatter, level=logging.INFO)
    logging.info(f'Starting {os.path.basename(__file__)}')

    logging.info(f'Loading IMGT alignments: {args.imgt_alignments} reference HLA bed: {args.ref_hla_bed},' + 
     f'reference HLA genotypes: {args.ref_hla_genotypes}')
    ref_bed_df = pd.read_csv(args.ref_hla_bed, header=None, sep='\t', 
        names=['chrom', 'start', 'end', 'score', 'strand'], index_col=3)
    ref_gt_df = pd.read_csv(args.ref_hla_genotypes, header=None, sep='\t', names=['genotype'], index_col=0)
    imgt_aln_gen = pickle.load(open(args.imgt_alignments, 'rb'))

    if os.path.isfile(args.alleles[0]):
        alleles_l = [x.strip() for x in open(args.alleles[0]).readlines()]
    else:
        alleles_l = args.alleles

    logging.info(f'Finding nearest IMGT alleles for {",".join(alleles_l)}')
    alleles_nearest = hla_nearest(alleles_l, imgt_aln_gen)
    df = aln2snp(imgt_aln_gen, ref_bed_df, ref_gt_df, alleles_nearest)

    logging.info(f'Writing {args.outfile}')
    df.to_csv(args.outfile, sep='\t', index=False)

if __name__ == "__main__":
    main()
