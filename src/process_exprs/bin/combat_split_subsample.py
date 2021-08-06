""" COMBAT batch correction subsampling stress test
Subsamples PACT samples are re-runs COMBAT.
"""
import sys
import pandas as pd
import patsy
import numpy as np
from datetime import datetime

sys.path.append('/home/csmith/git/combat.py/')
from combat2 import combat

sys.path.append('/home/csmith/git/bioinfo-tme/bin')
from combat_split import combat_split

if __name__ == "__main__":
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    print(dt_string)	

    # Config
    sites = ['ColonRectum', 'Breast', 'Lung', 'Skin', 'Ovary']

    # Set True to output hgnc symbols as rownames
    hgnc = True

    meta_f = '/media/nfs/data/References/xena/meta_pact_xena.tsv'
    meta_pact_f = '/home/csmith/csmith/tme/metadata_research_clinical_20210408.tsv'
    #exprs_xena_f = '/media/nfs/data/References/xena/subsets/TcgaTargetGtex_rsem_gene_tpm_COADREADOVLUADLUSCSKCMBRCA_1000.gz'
    exprs_xena_f = '/media/nfs/data/References/xena/TcgaTargetGtex_rsem_gene_tpm.gz'
    #exprs_xena_f = '/media/nfs/data/References/xena/subsets/TcgaTargetGtex_rsem_gene_tpm_SKCM.gz'
    exprs_pact_f = '/home/csmith/git/bioinfo-tme/pact/pact_rsem_tpm.exprs.tsv.gz'
    blacklist = [line.strip() for line in open("/home/csmith/git/bioinfo-tme/data/blacklist.txt", 'r')]
    hgnc_f = '/media/nfs/data/References/xena/hgnc_pact_xena.tsv'

    #exprs_out = '/media/nfs/data/References/xena/subsets/PactTcgaTargetGtex_rsem_gene_tpm'

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
    #print(exprs.shape)

    meta_pact = meta[meta._study == 'PACT']
    exprs_pact = pd.read_csv(exprs_pact_f, sep='\t', index_col=0)
    exprs_pact = exprs_pact.drop(blacklist, axis=1)

    from math import floor
    import random
    import numpy as np
    df = pd.DataFrame()
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

    # HGNC
    if hgnc:
        hgnc_df = pd.read_table(hgnc_f, index_col=0)
        hgnc_df = hgnc_df.drop(columns='_merge')
        exprs_m = exprs_m.rpow(2) 
        exprs_m = exprs_m.merge(hgnc_df, left_index=True, right_index=True)
        exprs_m['gene'] = np.where(exprs_m.gene.isna(), exprs_m.index, exprs_m.gene)
        exprs_m = np.log2(exprs_m.groupby('gene').sum())

    print('\nMerged')

    ## Metadata
    meta_m = pd.concat([meta_x, meta_pact])
    #sids = exprs.columns.to_list() + exprs_pact.columns.to_list()
    #meta_m = meta.loc[sids]
    print(meta_m.groupby(['_study', 'tumor_normal']).primary_site.value_counts())

    exprs_out = f'/media/nfs/data/Workspace/Users/csmith/tme/combat_stress'
    #for jj in np.linspace(0,1,11):
    # Combat
    print('\nCombat')
    df = pd.DataFrame()
    for ii in sites:
        meta_site = meta_m[meta_m.primary_site == ii]

        # Samples to drop
        tmp_pact = meta_site[meta_site._study == 'PACT']
        n_pact = tmp_pact.shape[0]
        n_samples = {floor(x) for x in np.linspace(0,1,11) * n_pact}
        # n = floor(tmp_pact.shape[0] * jj/100)
        for n in n_samples:
            p_kept = round((1 - n/n_pact) * 100)
            print(f'{ii}: {n_pact} ({p_kept}%) samples kept')

            tmp_drop = tmp_pact.sample(n)
            meta_site_n = meta_site.drop(tmp_drop.index) 
            exprs_site_n = exprs_m.loc[:, meta_site_n.index]
            
            edat_combat = combat_split(exprs_site_n, meta_site_n, batch_name = '_study', split_name = 'primary_site')
            if not edat_combat.empty:
                success = True
                outfile = f'{exprs_out}/exprs_combat_{ii}_hgnc_{p_kept}.parquet.gz'
                print(f'Writing {outfile}\n')
                exprs_site_n.to_parquet(outfile, compression='gzip')
            else:
                success = False
            n_kept = n_pact - n
            row = {
                'site': ii,
                'n_kept': n_kept,
                'p_kept': p_kept,
                'n_pact': n_pact,
                'ngenes': edat_combat.shape[0],
                'success': success,
                'fn': outfile,
            }
            tmp = pd.DataFrame([row])
            df = pd.concat([df, tmp], ignore_index=True)
    outfile = f'{exprs_out}/combat_stress.tsv'
    print(f'DONE: Writing {outfile}')
    df.to_csv(outfile, sep='\t', index=False)
        

