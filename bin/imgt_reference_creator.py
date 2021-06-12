#!/usr/bin/env python
''' 
IMGT alignment reference creator
'''
import sys
import pandas as pd
import pickle
import argparse
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

def parse_args(args=None):
    parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=__doc__)
    parser.add_argument('-i', '--imgt_path', required=True, 
        help='Path to IMGTHLA git repository (https://github.com/ANHIG/IMGTHLA)')
    parser.add_argument('-g', '--ref_hla_genotypes', required=True, help='Reference HLA genotypes')
    parser.add_argument('-b', '--ref_hla_bed', required=True, help='Reference HLA bed')
    parser.add_argument('-o', '--outfile', default='imgt_aln_gen.p', help='Output filename')
    
    return parser.parse_args(args)

def imgt_align(imgt_path, ref_bed_df, ref_gt_df):
    imgt_aln_gen = dict()
    for ii in ref_bed_df.index:
        ref_aln = None
        ref_gt = ref_gt_df.loc[ii].genotype
        fn = f'{imgt_path}/msf/{ii}_gen.msf'
        aln = AlignIO.read(fn, 'msf')
        for rec in aln:
            if rec.id == ref_gt:
                print(f'{ii} reference allele {ref_gt} found')
                ref_aln = MultipleSeqAlignment([rec])
        if not ref_aln:
            print(f'{ii} reference allele not found. Skipping...')
            continue
        for rec in aln:
            ref_aln.append(rec)
        imgt_aln_gen[ii] = ref_aln
    return imgt_aln_gen

def main(args=None):
    args = parse_args(args)
    ref_bed_df = pd.read_csv(args.ref_hla_bed, header=None, sep='\t',
        names=['chrom', 'start', 'end', 'score', 'strand'], index_col=3)
    ref_gt_df = pd.read_csv(args.ref_hla_genotypes, header=None, sep='\t', names=['genotype'], index_col=0)

    print('Generating IMGT alignments')
    imgt_aln_gen = imgt_align(args.imgt_path, ref_bed_df, ref_gt_df)

    print(f'Pickling to {args.outfile}')
    pickle.dump(imgt_aln_gen, open(args.outfile, 'wb'))

if __name__ == "__main__":
    main()