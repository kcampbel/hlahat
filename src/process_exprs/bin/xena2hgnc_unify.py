#!/usr/bin/env python
""" Generate Xena HGNC emat """
import pandas as pd
import logging
import os

formatter = '%(asctime)s:%(levelname)s:%(name)s:%(funcName)s: %(message)s'
logging.basicConfig(format=formatter, level=logging.INFO)
logging.info(f'Starting {os.path.basename(__file__)}')
 
infile = '/media/nfs/data/References/xena/TcgaTargetGtex_rsem_gene_tpm.gz'
outfile = '/media/nfs/data/References/xena/TcgaTargetGtex_rsem_hgnc_tpm.parquet.gz'
hgnc_f = '/home/csmith/git/bioinfo-fio/process_exprs/data/hgnc_pact_xena.tsv'

if os.path.exists(outfile):
    raise Exception(f'{outfile} exists. Aborting...')

hgnc = pd.read_csv(hgnc_f,sep='\t', index_col='ensembl_gene_id')
hgnc = hgnc.drop(columns='_merge')

logging.info(f'Reading {infile} and summarizing to gene from {hgnc_f}')
df_o = pd.read_csv(infile,sep='\t', index_col=0)
df = df_o
df.index = df.index.str.rsplit('.').str[0]
df = df.rpow(2) - 0.001
dfm = df.merge(hgnc, how='left', left_index=True, right_index=True)\
    .groupby('gene').apply(sum, axis=0)
dfm.index.name = 'hgnc_symbol'
dfm = dfm.drop(columns='gene')

logging.info(f'Writing {outfile}')
dfm.to_parquet(outfile, compression='gzip')
logging.info(f'{os.path.basename(__file__)} finished')
