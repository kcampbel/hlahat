import pandas as pd
import numpy as np

def seqz_join(seqz, other):
    """ Intersects a Sequenza seqz file with a pyranges object
    Args:
        seqz(pandas): Sequenza small.seqz.gz dataframe
        other(pr): pyranges object to join
        
    Returns:
        pyranges object
    """
    
    df = seqz[(seqz.chromosome == 'chr6') & (seqz.Bf > 0)]
    df = df.rename(columns={'chromosome': 'Chromosome'})\
      .assign(Start = df.position, End=df.position)
    pr_seqz = pr.PyRanges(df)
    pr_seqz_join = pr_seqz.join(other)
    
    out = pr_seqz_join.as_df()\
      .rename(columns = {
          'Chromosome': 'chr',
          'Start': 'start',
          'End': 'end'})
    return out 

def bcf_to_df(df, min_reads):
    df = df[~df.alt.str.contains('<*>')]
    #df['position'] = df.position.astype('int')
    nAD = df.nAD.str.split(',', expand=True).astype('int')
    tAD = df.tAD.str.split(',', expand=True).astype('int')
    df = df.assign(nREF = nAD[0], 
                   nALT = nAD[1], 
                   tREF = tAD[0],
                   tALT = tAD[1],
                   refCounts = tAD[0], 
                   varCounts = tAD[1])
    df = df.assign(
                   start = df.position - 1,
                   end = df.position,
                   normalCounts = df.nALT + df.nREF,
                   tumorCounts = df.varCounts + df.refCounts,
                   obsVaf = df.varCounts/(df.varCounts + df.refCounts),
                   nVaf = df.nALT/(df.nALT + df.nREF),
                  )\
      .drop(columns=['nAD', 'tAD', 'nGT', 'tGT'])#, 'nREF', 'nALT'])
    out = df[
        (
            df.normalCounts > min_reads
        )
    ]
    return(out)

def vaf_normalize_to_normal(df):
    """ Normalizes tumor VAF to normal VAF
    Normalizes tumor ref and alt ref counts to the normal for each SNP

    Args:
        df(pandas): requires tumor read counts tALT, tREF and normal read counts nALT, nREF

    Returns:
        pandas
    """
    tALT_norm, tREF_norm = df.tALT/df.nALT, df.tREF/df.nREF
    out = df.assign(tALT_norm=tALT_norm, tREF_norm=tREF_norm, obsVafNorm = tALT_norm/(tALT_norm + tREF_norm))
    return(out)

def flip_snps(df, vaf_col:str, prop:float = 0.5, new_col:str = 'obsVafNorm_flip'):
    """ SNP VAF flipper
    Flips VAF in a proportion of SNPs

    Args:
        df(pandas): pandas with columns: gene
        vaf_col(str): name of VAF column
        prop(float): proportion of SNPs to flip
        new_col(str): name of column with flipped VAFs
    
    Returns:
        pandas with new columns "vaf_col+'_flip'"(float) and "snp_flip"(bool)
    """
    samples_by_allele = (df.allele.value_counts()*prop).round(0).astype(int)

    out = df.copy()
    for allele in samples_by_allele.index:
        tmp = out[out.allele == allele]
        flipme = tmp.sample(samples_by_allele.loc[allele])
        flipme[new_col] = 1 - tmp[vaf_col]
        flipme['snp_flip'] = True
        out = out.drop(flipme.index, axis='index')
        out = pd.concat([out, flipme])
    out['snp_flip'] = out.snp_flip.fillna(False)
    out[new_col] = np.where(out[new_col].isna(), out[vaf_col], out[new_col])
    out.sort_index(inplace=True)
    return out



        

        

