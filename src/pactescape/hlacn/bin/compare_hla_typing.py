#!/usr/bin/env python
''' Compare HLA types from EPIC and hisat2-genotype
'''
import sys
import pandas as pd
import argparse
import re
import os
import csv
from collections import defaultdict, namedtuple, Counter
from pactescape.hlacn.imgt import hla_nearest

def parse_args(args=None):
    parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='EPIC-Hisat2 HLA type compare tool')
    parser.add_argument('query', help='optitype tsv')
    parser.add_argument('reference', help='hisat2-genotype text file')
    parser.add_argument('-o', '--outfile', default=sys.stdout)
    
    return parser.parse_args(args)

def hla_compare(query:dict, reference:dict):
    """ Epic to hisat hla compare tool 
    Compares a list of hla alleles to a reference, searching for substring matches. For 
    example, if query = A*02:01 and the reference = [A*02:01:01, A*03:01:01], A*02:01:01
    will be returned. If a match is not found, no alleles for that gene are returned. 

    Args:
        query(list): list of alleles
        reference(dict): {gene: [alleles]}
   
    Returns:
        List of alleles from query that match the reference
   """
    match = hla_nearest(query, reference)
    if len(match) == len(query):
        out = query
    else:
        # Number of alleles per locus in match
        match_n = Counter()
        for ii in match:
            gene = ii.split('*')[0]
            match_n[gene] += 1
        
        # Remove genes that don't have the same number of alleles as epic
        out = list()
        for k,v in match_n.items():
            epic_n = len([x for x in query if k in x])
            if v == epic_n:
                alleles = [x for x in match if k in x]
                out.extend(alleles)
    return(out)

def main():
    args = parse_args()
    #args = parse_args(['hlacn/test/hla-final.tsv', 'hlacn/test/hisatgenotype_normal_dna.report.tsv'])
    # Expected: A*02:01:01, A*03:01:01:05

    epic_df = pd.read_table(args.query,
        usecols = ['A1','A2','B1','B2','C1','C2'])
    epic_alleles = epic_df.iloc[0]

    hisat_d = defaultdict(list)
    al = namedtuple('Alleles', 'name')
    with open(args.reference) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            if int(row['rank']) > 2: #or not re.search('A|B|C', row['gene']):
                continue
            tup = al(row['allele'])
            hisat_d[row['gene']].append(tup)
            #gene = line.split('*')[0]
            #allele = line.strip()

    # Unit testing
#    epic_alleles = ['A*01:01', 'A*03:99', 'B*38:01', 'B*57:01']
#    hisat_d = {
#        'A': [al('A*02:01:01'), al('A*99:99:99')],
#        'B': [al('B*38:01:01'), al('B*57:01:01')],
#        'C': [al('C*01:01:01'), al('C*02:01:01')]
#    }
    out = pd.Series(hla_compare(epic_alleles, hisat_d))
    out.to_csv(args.outfile, index=False, header=False)

if __name__ == "__main__":
    main()
