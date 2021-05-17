#!/usr/bin/env python
""" Long format to expression matrix converter
Convert tsvs to a matrix with gene/transcript/symbol ids in rows and sample names in columns.
File name must be {sample_name}_tumor to parse the sample name for the column.
e.g. for RNA_PACT999_T_9999_tumor_rna.tar.gz, PACT999_T_9999 will be used as the column name.
""" 
import os
import sys
import pandas as pd 
import argparse

def main(counts:str, metric:str, fof:str, output:str, hgnc:bool):
    if fof:
        assert os.path.isfile(counts[0]), f'Counts must be a file of files if --fof is set.'
        reader = open(counts[0]).readlines()
        counts_l = [x.strip('\n') for x in reader]
    else:
        counts_l = counts

    df = pd.DataFrame()
    for ii in counts_l:
        sampleId = os.path.basename(ii).split('_tumor')[0].strip('RNA_')
        dat = pd.read_csv(ii, sep='\t')
        if hgnc:
            id_col = 'hgnc_symbol'
        else:
            id_col = dat.iloc[:,0].name
        tmp = dat.groupby(id_col)[metric].sum()
        tmp.name = sampleId
        tmp = tmp[tmp>0]
        df = df.merge(tmp, how='outer', left_index=True, right_index=True).fillna(0)
    df.to_csv(output, sep='\t')
    print(f'{len(counts_l)} samples processed')

if __name__ == '__main__':
    PARSER = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    PARSER.add_argument('metric', choices=['FPKM', 'TPM', 'counts'], help='Count metric')
    PARSER.add_argument('counts', nargs = '+', help='*genes.tsv or *isoforms.tsv files to process. File of file accepted with -f flag.')
    PARSER.add_argument('-f', '--fof', action='store_true', help='file of files to process')
    PARSER.add_argument('-o', '--output', default=sys.stdout, help='Output file')
    PARSER.add_argument('--hgnc', action='store_true', help='Return HGNC symbol as gene id')
    
    if len(sys.argv) == 1:
        PARSER.print_help()
        sys.exit()

    ARGS = PARSER.parse_args()
    #ARGS = PARSER.parse_args(['FPKM', 
    # 'pact/RNA_PACT004_T_560351F_tumor_rna.genes.tsv', 
    # 'pact/RNA_PACT384_T_021F_tumor_rna.genes.tsv'
    # '--hgnc'
    #])

    main(**vars(ARGS))
