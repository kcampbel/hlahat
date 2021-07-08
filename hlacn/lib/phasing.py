from io import StringIO
import pandas as pd
import numpy as np
from scipy.stats import chisquare, chi2_contingency

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

def bcftools_cmd(nBam, tBam, posBed, genome_fa, posTsv):
    posFormat = "'%CHROM\t%POS\t%REF\t%ALT[\t%GT][\t%AD]\n'" 
    cmd = ["bcftools", "mpileup",
            "-d", "5000", "-xBQ0", "--ignore-RG",
            "-a", "FORMAT/AD",
            "-R", posBed,
            "-f", genome_fa,
            nBam, tBam, "|",
            "bcftools", "norm", "-m-", "|",
            "bcftools", "query", "-uHf", posFormat, ">",
            posTsv]
    cmd = ' '.join(cmd)
    return(cmd)

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

def fitSequenza_cmd(specimen_id:str, varsFile:str, outputDir:str, altMode:str, sequenzaModelRData:str, sequenzaTools:str):
    cmd = [
        'fitSequenzaModel_hla.R',
        '-i', specimen_id,
        '-v', varsFile,
        '-o', outputDir,
        '-a', altMode,
        '-r', sequenzaModelRData,
        '-t', sequenzaTools
    ]
    return(cmd)

def chisq2dict(obs, CNt:int = None):
    """ Chisquare results to dictionary
    Args:
        obs(pandas): observed frequencies with alleles in rows and Mt in columns
        
    Returns:
        dict of {'c' = chisquare, 'p': p-value} 
        
    """
    if len(obs.shape) > 1:
        result = chi2_contingency(obs)
        out = {
            'c': result[0],
            'p': result[1],
            #'dof': result[2],
            #'expected': result[3]
            }
    else:
        
        result = chisquare(obs) 
        out = {
            'c': result[0],
            'p': result[1],
            }
        
    return out

def vaf_normalize_to_normal(df, method:str = 'snp_count_ratio'):
    out = pd.DataFrame()
    if method == 'snp_count_ratio':
        out = df
        tALT_norm, tREF_norm = df.tALT/df.nALT, df.tREF/df.nREF
        out = out.assign(tALT_norm=tALT_norm, tREF_norm=tREF_norm, obsVafNorm = tALT_norm/(tALT_norm + tREF_norm))
    else:
        df_grp = df.groupby(['Start_b', 'End_b'])
        if method == 'count_sum':
            tmp = df_grp.agg(nREF_win=('nREF', 'sum'), nALT_win=('nALT', sum), snps_per_window=('nALT', 'count')).reset_index()
            tmp['nVafWindow'] = tmp.nALT_win/(tmp.nALT_win + tmp.nREF_win)
        elif method == 'nVaf_median':
            tmp = df_grp.agg(nVafWindow = ('nVaf', 'median'), snps_per_window=('nALT', 'count')).reset_index()
        elif method == 'nVaf_mean':
            tmp = df_grp.agg(nVafWindow = ('nVaf', 'mean'), snps_per_window=('nALT', 'count')).reset_index()
        elif method == 'nVaf_weighted_mean':
            tmp = df_grp.apply(lambda x: np.average(x.nVaf, weights=x.normalCounts)).reset_index(name='nVafWindow')
        else:
            print(f'{method} not found')
            return None
        out = df.merge(tmp)
        obsVafNorm = out.obsVaf + (0.5 - out.nVafWindow)
        obsVafNorm = np.where(obsVafNorm > 1, 1, obsVafNorm)
        obsVafNorm = np.where(obsVafNorm < 0, 0, obsVafNorm)
        out['obsVafNorm'] = obsVafNorm
    return(out)

