""" COMBAT batch correction
"""
import sys
import pandas as pd
sys.path.append('/home/csmith/git/combat.py/')
from combat2 import combat
import patsy
import numpy as np
from datetime import datetime

def combat_split(exprs, meta, batch_name, split_name, todict:bool = False):
    if todict:
        out = dict()
    else:
        out = pd.DataFrame()
    for ii in meta[split_name].unique():
        print(ii)
        mf = meta[meta[split_name] == ii]
        batch = mf[batch_name]
        samples = mf.index
        tmp = exprs[samples]
        #tmp = tmp.loc[tmp.apply(lambda x: all(x > np.log2(0 + 0.001)), axis=1)]

        if mf.study_tumor_normal.str.contains('TCGANormal').any():
            formula = "~ tumor_normal"
        else:
            formula = '~ 1'
        print(formula)
        mod = patsy.dmatrix(formula, mf, return_type="dataframe")
        edat = combat(tmp, batch, mod)
        
        if (edat.isna().apply(sum) != 0).any():
            print(f'WARNING: {ii} failed COMBAT. Skipping...')
            continue
        if todict:
            out[ii] = edat
        else:
            out = pd.concat([out, edat], axis=1)
    return out

now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
print(dt_string)	

# Config
#sites = ['ColonRectum', 'Lung', 'Skin', 'Breast', 'Ovary']
sites = ['Breast']
meta_f = '/media/nfs/data/References/xena/meta_pact_xena.tsv'
#TcgaTargetGTEX_phenotype_utf8.txt'
meta_pact_f = '/home/csmith/csmith/tme/metadata_research_clinical_20210408.tsv'

#exprs_xena_f = '/media/nfs/data/References/xena/subsets/TcgaTargetGtex_rsem_gene_tpm_COADREADOVLUADLUSCSKCMBRCA_1000.gz'
exprs_xena_f = '/media/nfs/data/References/xena/TcgaTargetGtex_rsem_gene_tpm.gz'
#exprs_xena_f = '/media/nfs/data/References/xena/subsets/TcgaTargetGtex_rsem_gene_tpm_SKCM.gz'
exprs_pact_f = '/home/csmith/git/bioinfo-tme/pact/pact_rsem_tpm.exprs.tsv.gz'
blacklist = [line.strip() for line in open("/home/csmith/git/bioinfo-tme/pact/blacklist.txt", 'r')]

exprs_out = '/media/nfs/data/References/xena/subsets/PactTcgaTargetGtex_rsem_gene_tpm'

# Metadata
meta = pd.read_csv(meta_f, sep='\t', index_col=0)
meta = meta[
    (meta.primary_site.isin(sites)) &
    (meta._sample_type != 'Cell Line')
]
print(f'Sites: {",".join(sites)}')
print('\nXena')
print(meta.groupby(['_study', 'tumor_normal']).primary_site.value_counts())
print(f'{meta.shape[0]} samples total')

# Expression matrices
meta_x = meta[meta._study != 'PACT']
sids = meta_x.index.to_list()
exprs = pd.read_csv(exprs_xena_f, sep='\t', index_col=0, usecols=['sample'] + sids)
print(exprs.shape)

#sids = meta[meta._study == 'PACT'].index.to_list()
#meta_pact = pd.read_table(meta_pact_f, index_col=0)
meta_pact = meta[meta._study == 'PACT']
exprs_pact = pd.read_csv(exprs_pact_f, sep='\t', index_col=0)
exprs_pact = exprs_pact.drop(blacklist, axis=1)

# Harmonize
exprs_pact = exprs_pact.loc[:, exprs_pact.columns.isin(meta.index)]
meta_pact = meta_pact.loc[exprs_pact.columns]
exprs_pact = np.log2(exprs_pact + 0.001)

print('\nPACT')
print(f'{meta_pact.shape[0]} samples total')

# Merge
## expression
# Remove trailing .X from gene Ids
exprs.index = exprs.index.str.rsplit('.').str[0]
exprs_m = pd.concat([exprs, exprs_pact], axis=1).fillna(np.log2(0 + 0.001))

print('\nMerged')
print(exprs_m.shape)

## Metadata
meta_m = pd.concat([meta_x, meta_pact])
#sids = exprs.columns.to_list() + exprs_pact.columns.to_list()
#meta_m = meta.loc[sids]
print(meta_m.groupby(['_study', 'tumor_normal']).primary_site.value_counts())

# Combat
print('\nCombat')
for ii in sites:
    meta_site = meta_m[meta_m.primary_site == ii]
    exprs_site = exprs_m.loc[:, meta_site.index]
    edat_combat = combat_split(exprs_site, meta_site, batch_name = '_study', split_name = 'primary_site')
    outfile = f'{exprs_out}/exprs_combat_{ii}.tsv.gz'
    print(f'Writing to {outfile}')
    edat_combat.to_csv(outfile, sep='\t')

now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
print(dt_string)	


