#!/usr/bin/env python
''' Hisat-genotype report to tsv
'''
import sys
import pandas as pd
import argparse
import re
import os

def parse_args(args=None):
    parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='Hisat-genotype report to tsv conversion tool')
    parser.add_argument('infile')
    parser.add_argument('-o', '--outfile', default=sys.stdout)
    
    return parser.parse_args(args)

def main():
    args = parse_args()
#    args = parse_args(
#        [
#            'test/PACT999_T_999_normal_dna.report'
#        ]
#    )
    df = pd.DataFrame()
    with open(args.infile, 'rt') as fh:
        for line in fh:
            if 'pairs are aligned' in line:
                rank = 0
                tmp = line.lstrip().split(' ')
                reads, pairs = tmp[0], tmp[3]
            if 'ranked' in line:
                rank += 1
                tmp = line.lstrip().split(' ')
                allele = tmp[2]
                abundance = tmp[4].split('%)')[0]
                row = {
                    'allele': allele,
                    'gene': allele.split('*')[0],
                    'reads': reads,
                    'pairs': pairs,
                    'abundance': abundance,
                    'rank': rank,
                }
                df = pd.concat([df, pd.DataFrame([row])])
    df = df.sort_values(['gene', 'rank'])
    df.to_csv(args.outfile, sep='\t', index=False)

if __name__ == '__main__':
    main() 




       

