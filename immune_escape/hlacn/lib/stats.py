import pandas as pd
import numpy as np
import logging
from scipy.stats import chisquare, chi2_contingency
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm
from scipy.stats import mannwhitneyu

def make_stats_df(df, hla_types:dict, vaf_col:str = 'obsVafNorm'):
    """ SNP AF stats input dataframe creator
    Generates a df with 'success' and 'fail' counts for binomial GLMs, and vaf column for non-parametric stats.
    The boolean column 'snp_flip' is used to determine whether ALT,REF counts should be assigned to 'success' or 'fail'.
    HLA typing is required to label cases where one allele matches the reference and thus contains no SNPs.

    Args:
        df(pandas): contains columns gene, snp_flip, tALT_norm, tREF_norm, allele, success, fail, vaf where rows are SNPs
        hla_types(dict): hla genotypes where {gene_name: [allele1, allele2]}
        vaf_col(str): name of VAF column for non-parametric statistics

    Returns:
        pandas with gene, allele, success, fail, vaf
    """
    out = pd.DataFrame()
    for Index, row in df.iterrows():
        gene = row.gene
        alleles = np.array(hla_types[gene])
        if row.snp_flip:
            success, fail = row.tREF_norm, row.tALT_norm
            group = alleles[alleles != row.allele][0]
            vaf = 1 - row[vaf_col]
        else:
            success, fail = row.tALT_norm, row.tREF_norm
            group = alleles[alleles == row.allele][0]
            vaf = row[vaf_col]
        rec = {
            'gene': gene,
            'allele': group,
            'success': success,
            'fail': fail,
            'vaf': vaf,
        }
        tmp = pd.DataFrame([rec])
        out = pd.concat([out, tmp])
    return(out)

def compare_vaf(df, min_snps_per_allele:int = 5):
    """ Compare AF between alleles
    Compares counts (binomial GLM) or vaf (Mann-Whitney U test) between alleles for each gene in a pandas df.
    See make_stats_df to generate input.

    Args:
        df(pandas): contains columns gene, allele, success, fail, vaf where rows are SNPs
        min_snps_per_allele(int): minimum number of snps to use nonparametric test
    
    Returns:
        pandas with gene, tstatistic, df, p-value
    """
    def no_stats(gene, test):
        row = {
            'gene': gene,
            'method': test,
            'tstat': None,
            'pvalue': None
        }
        return row
    
    out = pd.DataFrame()
    for gene in df.gene.unique():
        df_g = df[df.gene == gene]
        nsnps = df_g.groupby('gene').allele.value_counts()
        # Skip if homozygous
        if df_g.allele.unique().size < 2:
            logging.warning(f'{gene} is reported as homozygous {df_g.allele.unique()}. Skipping...')
            row = no_stats(gene, None)
            continue

        # Do parametric test if min_snps_per_allele is reached, otherwise non-parametric
        if (nsnps > min_snps_per_allele).any():
            test = 'binomialglm'
            logging.info(f'Comparing AF for {gene} with {test}')
            tstat, pvalue, dof = stats_vaf(df_g, test)
        else:
            test = 'mannwhitneyu'
            logging.info(f'{gene} has SNP counts below the {min_snps_per_allele} snp minimum\n{nsnps}')
            logging.info(f'Comparing AF for {gene} with {test}')
            tstat, pvalue, dof = stats_vaf(df_g, test)
        row = {
            'gene': gene,
            'method': test,
            'tstat': tstat,
            'df': dof,
            'pvalue': pvalue
        }
        tmp = pd.DataFrame([row])
        out = pd.concat([out, tmp])
    out = out.set_index('gene')
    return out

def stats_vaf(df, test:str, vaf_col:str = 'vaf'):
    """ Methods for comparing AF
    Compares counts (binomial GLM) or vaf (Mann-Whitney U test) between alleles 

    Args:
        df(pandas): contains columns allele, success, fail, vaf where rows are SNPs
        min_snps_per_allele(int): minimum number of snps to use nonparametric test
    
    Returns:
        tstat, pvalue, df
    """
    if test == 'mannwhitneyu':
        alleles = df.allele.unique()
        x = df[df.allele == alleles[0]][vaf_col]
        y = df[df.allele == alleles[1]][vaf_col]
        result = mannwhitneyu(x, y)
        tstat, pvalue, dof = result[0], result[1], None
    elif test == 'binomialglm':
        mod = smf.glm('success + fail ~ allele', family=sm.families.Binomial(), data=df).fit()
        A = np.identity(len(mod.params))
        A = A[1:,:]
        result = mod.f_test(A)
        tstat, pvalue, dof = result.statistic.flatten()[0], result.pvalue, (result.df_num, result.df_denom)
    return(tstat, pvalue, dof)
    
def get_major_allele(df, y, stat = 'mean'):
    tmp = df.groupby(['gene', 'allele'])[y].agg(stat).reset_index()
    max_vaf = tmp[y].max()
    major_allele = tmp[(tmp[y] == max_vaf)].allele.iloc[0]
    minor_allele = tmp[(tmp[y] != max_vaf)].allele.iloc[0]
    return major_allele, minor_allele

def allele_descr_stats(df, y, stat:str = 'mean'):
    """ Allele SNP descriptive stats
    Calculates descriptive stats of 'y' for each record by gene in the df and determines the allele 
    with the higher and lower statistic 'stat'. Useful for determining which allele
    has the higher copy number, vaf, etc.

    Args:
        df(pandas): contains columns gene, allele, y where rows are SNPs
        y(str): column name for calculating statistics
        stat(str): statistic to use (eg median, mean)

    Returns:
        pandas with gene, major allele, and minor allele, and descriptive stats
    """
    dstats = df.groupby(['gene', 'allele'])[y].describe()\
        .reset_index()\
        .rename(columns={'50%': 'median'})
    
    major_minor = pd.DataFrame()
    for gene in df.gene.unique():
        tmp = df[df.gene == gene]
        if tmp.allele.unique().size != 2:
            row = {
                'gene': gene,
                'major_allele': None,
                'minor_allele': None,
            }
            major_minor = pd.concat([major_minor, pd.DataFrame([row])], ignore_index=True)
        else:
            major_allele, minor_allele = get_major_allele(tmp, y, stat)
            row = {
                'gene': gene, 
                'major_allele': major_allele, 
                'minor_allele': minor_allele,
            }
                
            major_minor = pd.concat([major_minor, pd.DataFrame([row])], ignore_index=True)
    out = dstats.merge(major_minor)
    return out

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