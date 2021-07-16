import pandas as pd
import os
import logging

def counts2mat(counts:list, metric:str, gene_name:str, exprs, force:bool):
    if exprs.index.name != gene_name:
        raise Exception(f'exprs gene column name {exprs.index.name} does not match argument --gene_name {gene_name}. Aborting...')
    # Add new samples
    added, skipped = list(), list()
    for ii in counts:
        sampleId = os.path.basename(ii).split('_tumor')[0].strip('RNA_')
        if exprs.columns.isin([sampleId]).any():
            if force:
                logging.info(f'--force enabled. Dropping {sampleId} from exprs')
                exprs = exprs.drop(sampleId, axis=1)
            else:
                logging.warning(f'{sampleId} already exists in exprs. Set --force to overwrite.')
                skipped.append(ii)
                continue
        dat = pd.read_csv(ii, sep='\t')

        # Sum rows by gene name and remove genes with zero expression
        tmp = dat.groupby(gene_name)[metric].sum()
        tmp.name = sampleId
        tmp = tmp[tmp > 0]

        # Merge with exprs matrix
        logging.info(f'Adding {sampleId} to exprs with {exprs.shape[1]} samples')
        exprs = exprs.merge(tmp, how='outer', left_index=True, right_index=True).fillna(0)
        added.append(ii)

    logging.info(f'{metric} from {len(added)} samples added and {len(skipped)} skipped.')
    if added:
        return exprs
    else:
        return None

def manifest_to_meta(manifest:dict):
    """ Converts EPIC manifest fields to a dict
    """
    pi = manifest['pipeline']['patient_info']
    row = {
        'specimen_id': pi['patient'],
        'pact_patient_id': pi['patient'].split('_')[0],
        'patient_id': pi['patient.id'],
        'sample_name': pi['sample.id'],
        'date_of_birth': pi['dob'],
        'study_id': pi['study.id'],
        'tcga_study_code': pi['patient.tumorType'],
    }
    return(row)
