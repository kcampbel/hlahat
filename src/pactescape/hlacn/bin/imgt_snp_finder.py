#!/usr/bin/env python
""" IMGT SNP finder

Determines SNPs for input HLA alleles using a custom IMGT alignment where the reference
genome HLA types are the reference alleles. See imgt_reference_creator.py to produce the custom 
alignment.
"""
import sys
import pandas as pd
import pickle
import argparse
import os
import logging
from pactescape.hlacn.imgt import arr2snp, aln2snp, hla_nearest
from commonLib.lib.fileio import package_file_path
from pactescape.hlacn import data

def parse_args(args=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument('alleles', nargs='+', help='Query alleles or a tsv')
    parser.add_argument('-i', '--imgt_alignments', 
        default=package_file_path(data, 'imgt_aln_gen.p'),
        help='IMGT alignments pickle input')
    parser.add_argument('-g', '--ref_hla_genotypes', 
        default=package_file_path(data, 'hla_genotypes_hs37d5.tsv'),
        help='Reference HLA genotypes')
    parser.add_argument('-b', '--ref_hla_bed',
        default=package_file_path(data, 'hla_hs37d5.bed'), 
        help='Reference HLA bed')
    parser.add_argument('-o', '--outfile', default=sys.stdout)
    
    return parser.parse_args(args)


def main(args=None):
    args = parse_args(args)
#    args = parse_args(
#        [
#            '--imgt_alignments', '/home/csmith/git/bioinfo-fio/immune_escape/hlacn/src/hlacn/data/imgt_aln_gen.p',
#            '-g', '/home/csmith/git/bioinfo-fio/immune_escape/hlacn/src/hlacn/data/hla_genotypes_hs37d5.tsv',
#            '-b', '/home/csmith/git/bioinfo-fio/immune_escape/hlacn/src/hlacn/data/hla_hs37d5.bed',
#            '-o', '/tmp/snps.tsv',
#            'A*02:01',
#        ]
#    )
    logging.info(f'Starting {os.path.basename(__file__)}')

    logging.info(f'Loading IMGT alignments: {args.imgt_alignments} reference HLA bed: {args.ref_hla_bed},' + 
     f'reference HLA genotypes: {args.ref_hla_genotypes}')
    ref_bed_df = pd.read_csv(args.ref_hla_bed, header=None, sep='\t', 
        names=['chrom', 'start', 'end', 'score', 'strand'], index_col=3)
    ref_gt_df = pd.read_csv(args.ref_hla_genotypes, header=None, sep='\t', names=['genotype'], index_col=0)
    imgt_aln_gen = pickle.load(open(args.imgt_alignments, 'rb'))

    if os.path.isfile(args.alleles[0]):
        alleles_df = pd.read_csv(args.alleles[0], sep='\t')
        alleles_l = alleles_df[alleles_df.gene.isin(['A','B','C'])].alleles
    else:
        alleles_l = args.alleles

    logging.info(f'Finding nearest IMGT alleles for {",".join(alleles_l)}')
    alleles_nearest = hla_nearest(alleles_l, imgt_aln_gen)
    df = aln2snp(imgt_aln_gen, ref_bed_df, ref_gt_df, alleles_nearest)
    logging.info(f'Writing {args.outfile}')
    df.to_csv(args.outfile, sep='\t', index=False)

if __name__ == "__main__":
    main()
