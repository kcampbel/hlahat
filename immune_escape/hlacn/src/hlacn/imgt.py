from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
import pandas as pd
import numpy as np
import re
import os
import logging

def arr2snp(arr, alleles:list, gene_start:int = 0, chromosome:str = "chr6"):
    """ Numpy array SNP parser
    Parses SNPs and returns their position, adding gene_start if liftover is required.  
    The first row in arr is the reference. Insertions and deletions in the alignment are 
    not included in the output.

    .param arr(numpy): Numpy array with rows as sequences and columns as positions
    .param alleles(list): Names of alleles in rows of arr
    .param gene_start(int): Position in reference where gene start is located
    .param chromosome(str): Name of chromosome for reporting

    .rettype pandas with SNPs for each row in arr
    """
    # Remove insertions
    ins = np.nonzero(arr[0] != '-')[0]
    arr = arr[:, ins]
    
    # Find SNPs
    diffs = np.nonzero(arr != arr[0,:])
    
    out = pd.DataFrame()
    for ii in np.transpose(diffs):
        snp = arr[ii[0],ii[1]]
        if snp == '-':
            continue
        row = {
            'chromosome': chromosome,
            'position': ii[1] + gene_start,
            'snp': snp,
            'allele': alleles[ii[0]],
            'allele.ref': alleles[0],
            'pos_imgt': ii[1],
            'base.ref': arr[0, ii[1]]
        }
        tmp = pd.DataFrame([row])
        out = pd.concat([out, tmp])
    return out

def hla_nearest(query:list, imgt:dict):
    """ HLA allele nearest match
    Searches imgt[locus] = alleles for an exact match. If no match is found, search for a more
    specific match (e.g. query=A*02:01, IMGT match=A*02:01:01). If a more specific match is 
    not found, reduce the field resolution until a match is found 
    (e.g. query=A*02:01:99, IMGT match = A*02:01:01)
    
    Args:
        query(list): query alleles
        imgt(dict): {locus:[alleles]}
        
    Returns:
        List of nearest matches
    """
    out = list()
    for jj in query:
        locus = jj.split('*')[0]
        imgt_alleles = [x.name for x in imgt[locus]]
        # Search for exact match
        exact_match = [jj == x for x in imgt_alleles]
        if any(exact_match):
            result = imgt_alleles[exact_match.index(True)]
        else:
            # Find a match with more fields
            nearest_longer = [x for x in imgt_alleles if jj in x]
            if nearest_longer:
                result = nearest_longer[0]
                logging.info(f'{jj} -> {result} nearest match')
            else:
                # Reduce to two-field resolution until a match is found
                field_range = range(1, len(jj.split(':')) - 1)
                if len(field_range) == 0:
                    result = None
                else:
                    for nfields in field_range:
                        pattern = jj.rsplit(':', nfields)[0]
                        nearest_same = [x for x in imgt_alleles if pattern in x]
                        if any(nearest_same):
                            result = nearest_same[0]
                            logging.info(f'{jj} -> {result} nearest match')
                            break
                        else:
                            result = None
        if not result:
            logging.warning(f'{jj} not found in reference. Skipping...')
            continue
        else:
            out.append(result)
    return(out)

def aln2snp(imgt_aln_gen, ref_bed_df, ref_gt_df, alleles):
    """ Parse SNPs from IMGT alignment
    Args
        imgt_aln_gen(pickle): Alignments in {'locus':MultipleSequenceAlignment} format, where
            the first record in the MSA is the reference alignment
        ref_bed_df(pandas): Bed of HLA loci with strand information
        ref_gt_df(pandas): Reference genotypes
        alleles(list): Query alleles

    Returns:
        pandas of SNPs and HLA types
    """
    # Remove loci from imgt alignments not in query
    loci = {x[0] for x in alleles}
    for k in loci: 
        if k not in imgt_aln_gen.keys():
            del imgt_aln_gen[k]
            
    
    out = pd.DataFrame()
    for k,v in imgt_aln_gen.items():
        # Subset alignment if alleles is specified, otherwise generate SNPs for all records
        allele_ref = ref_gt_df.loc[k].genotype
        alleles_locus = [x for x in alleles if re.match(k, x)]
        alleles_sample = [x for x in alleles_locus if x != allele_ref]
        alleles_all = [allele_ref] + alleles_sample

        alns_index_d = dict()
        for Count, ii in enumerate(v):
            alns_index_d[ii.name] = Count
        alns_l = list()
        for ii in alleles_all:
            alns_index = alns_index_d[ii]
            alns_l.append(v[alns_index, :])
        alns = alns_l

        # Reverse complement if locus is on minus strand
        if ref_bed_df.loc[k].strand == '-':
            recs = list()
            for ii in alns:
                ii.seq = ii.seq.reverse_complement()
                recs.append(ii)
            alns = MultipleSeqAlignment(recs)

        # Generate array and find SNPs
        arr = np.array(alns, dtype='str')
        gene_start = ref_bed_df.loc[k].start
        df = arr2snp(arr, alleles=alleles_all, gene_start=gene_start)
        out = pd.concat([out, df], ignore_index=True)
    out['gene'] = out.allele.str.split('*').str[0]
    return(out)
