#!/usr/bin/env python
""" HLA SNP Analysis Tool 

Calls SNPs from Tumor/Normal bams and generates copy number calls for HLA alleles using a Sequenza
copy number model. See imgt_snp_finder.py to produce the input SNPs.
"""
import sys
import os
import pandas as pd
import argparse
import subprocess as sb
import logging
import numpy as np
from commonLib.lib.fileio import package_file_path
from commonLib.lib.munge import merge_indicator_rename
from pactescape.hlacn.command import bcftools_cmd, fitSequenza_cmd
from pactescape.hlacn.munge import bcf_to_df, vaf_normalize_to_normal, flip_snps
import pactescape.hlacn

def parse_args(args=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=__doc__)
    parser.add_argument('-i', '--specimen_id', required=True, help='Specimen id')
    parser.add_argument('-n', '--normal_dna_bam', required=True, help='Normal BAM')
    parser.add_argument('-t', '--tumor_dna_bam', required=True, help='Tumor BAM')
    parser.add_argument('-g', '--genome_fasta', required=True, help='Reference genome FASTA')
    parser.add_argument('-p', '--snp_tsv', required=True, help='IMGT SNPs tsv')
    parser.add_argument('-m', '--sequenzaModelRData', required=True, help='Sequenza model RData file')
    parser.add_argument('-s', '--solution', type=str, default=1, help='Sequenza solution to analyze')
    parser.add_argument('-o', '--outDir', required=True, help='Output dir')
    
    return parser.parse_args(args)

def snp_qc_munger(df, merge:bool= False, labels:list = None):
    if merge:
        cnt = df.groupby('gene')._merge.value_counts()
        cnt.name = 'snps'
        cnt.index.set_names(['gene', 'set'], inplace=True)
        cnt = pd.DataFrame(cnt).reset_index()
        cnt['set'] = merge_indicator_rename(cnt.set, labels)
        cnt = cnt.pivot_table(values='snps', index='set', columns='gene')
    else:
        cnt = df.gene.value_counts()
        cnt.index.rename('gene', inplace=True)
        cnt.rename('snps', inplace=True)
        cnt = pd.DataFrame(cnt)
    return(cnt)

def main():
    args = parse_args()
#    args = parse_args(
#    [
#        '-i', 'PACT056_T_196454',
#        '-n', '/media/nfs/data/Workspace/Users/csmith/nf-GATK_Exome_Preprocess/asterand/PACT056_N_196450/PACT056_N_196450/PACT056_N_196450_normal_dna-final.bam',
#        '-t' ,'/media/nfs/data/Workspace/Users/csmith/nf-GATK_Exome_Preprocess/asterand/PACT056_T_196454/PACT056_T_196454/PACT056_T_196454_tumor_dna-final.bam',
#        '-g', '/media/nfs/data/References/hs37d5.chr/hs37d5.chr.fa',
#        '-p', '/home/csmith/git/bioinfo-fio/immune_escape/hlacn/test_data/PACT056_T_196454/PACT056_N_196450_snps.tsv',
#        '-m', '/home/csmith/git/bioinfo-fio/immune_escape/hlacn/test_data/PACT056_T_196454/sequenzaModel.RData',
#        '-o', '/tmp/snp2allele',
#        #'-s', '2'
#    ])
    logging.info(f'Starting {os.path.basename(__file__)}')
    logging.info(args)

    specimen_id = args.specimen_id
    nBam = args.normal_dna_bam
    tBam = args.tumor_dna_bam
    snp_f = args.snp_tsv
    genome_fasta = args.genome_fasta
    sequenzaModelRData = args.sequenzaModelRData
    outDir = f'{args.outDir}'
    logDir = f'{outDir}/log'

    missing = list()
    for ii in [tBam, nBam, snp_f, sequenzaModelRData]:
        if not os.path.exists(ii):
            missing.append(ii)
    if missing:
        raise FileNotFoundError(missing)

    if not os.path.exists(logDir):
        os.makedirs(logDir, exist_ok=True)

    snp_qc, step = pd.DataFrame(), 0
    # BCFTOOLS
    ## Make SNP bed
    snp_df = pd.read_table(snp_f)
    cnt = snp_qc_munger(snp_df)
    cnt['set'] = 'IMGT SNPs'
    cnt = cnt.pivot_table(values='snps', index='set', columns='gene')
    logging.info(f'IMGT SNPs: \n{cnt}')
    cnt['step'] = step
    snp_qc = pd.concat([snp_qc, cnt])

    snp_df = snp_df.rename(columns={'snp': 'alt', 'base.ref': 'ref'})
    snp_bed = snp_df.assign(start=snp_df.position-1).loc[:, ['chromosome', 'start', 'position']]
    snp_bed_f = f'{outDir}/{specimen_id}_snps.bed'
    snp_bed.to_csv(snp_bed_f, sep='\t', index=False, header=False)

    ## Run BCFtools
    posTsv = f'{outDir}/{specimen_id}_pos.tsv'
    cmd = bcftools_cmd(nBam, tBam, snp_bed_f, genome_fasta, posTsv)
    logging.info(cmd)
    logfile = f'{logDir}/bcftools.log'
    handle = open(logfile, 'wt')
    job = sb.run(cmd, stdout=handle, stderr=sb.STDOUT, shell=True, check=True)
    handle.close()

    # Filter SNPs
    ## Intersect with IMGT SNPs
    pos_df = pd.read_table(posTsv, skiprows=1, names=['chromosome', 'position', 'ref', 'alt', 'nGT', 'tGT', 'nAD', 'tAD'])

    tmp = bcf_to_df(pos_df, min_reads=0)
    pos_m = tmp.merge(snp_df, how='outer', indicator=True).drop_duplicates()
    labels = ['bcftools only', 'IMGT only', 'bcftools and IMGT']
    cnt = snp_qc_munger(pos_m, merge=True, labels=labels)
    logging.info(f'Merging BCFtools and IMGT SNPs\n{cnt}')

    ### Munge for QC df
    cnt = cnt.loc['bcftools and IMGT']
    cnt.rename('snps', inplace=True)
    cnt = pd.DataFrame(cnt)
    cnt['set'] = 'bcftools and IMGT'
    cnt = cnt.pivot_table(values='snps', index='set', columns='gene')
    step = step + 1
    cnt['step'] = step
    snp_qc = pd.concat([snp_qc, cnt])

    ## Minimum normal DP
    min_normal = 20
    pos_imgt = pos_m[pos_m._merge == 'both']
    nsnps = pos_imgt.shape[0]
    pos_imgt = pos_imgt[pos_imgt.normalCounts > min_normal]
    cnt = snp_qc_munger(pos_imgt)
    logging.info(f'Removing {nsnps - pos_imgt.shape[0]} SNPs with less than {min_normal} reads in normal')
    cnt['set'] = 'Normal read min'
    cnt = cnt.pivot_table(values='snps', index='set', columns='gene')
    step = step + 1
    cnt['step'] = step
    snp_qc = pd.concat([snp_qc, cnt])
    
    ## Drop shared SNPs between allele 1 and allele 2
    snp_uniq = (pos_imgt.groupby(['position']).pos_imgt.count() == 1).reset_index()
    snp_uniq = snp_uniq[snp_uniq.pos_imgt].drop('pos_imgt', axis=1)
    nsnps = pos_imgt.shape[0]
    pos_imgt = pos_imgt.merge(snp_uniq)
    pos_imgt = pos_imgt.astype({'start': 'int32', 'end': 'int32'})

    cnt = snp_qc_munger(pos_imgt)
    cnt['set'] = 'Shared SNPs'
    cnt = cnt.pivot_table(values='snps', index='set', columns='gene')
    logging.info(f'Dropping {nsnps - pos_imgt.shape[0]} SNPs shared between alleles')
    step = step + 1
    cnt['step'] = step
    snp_qc = pd.concat([snp_qc, cnt])

    # Normalize tumor VAF to normal
    ## Drop SNPs with no REF counts
    tmp = pos_imgt[
        (pos_imgt.nREF == 0) |
        (pos_imgt.tREF == 0)
    ]
    if not tmp.empty:
        nempty = tmp.shape[0]
        cols = ["allele", "chromosome", "position", "obsVaf", "nREF", "tREF", "tALT", "nALT"]
        logging.warning(f'{nempty} SNP(s) dropped due to zero ref counts: \n{tmp[cols]}')
        pos_imgt = pos_imgt.drop(tmp.index)
    cnt = snp_qc_munger(pos_imgt)
    cnt['set'] = 'Zero REF counts'
    cnt = cnt.pivot_table(values='snps', index='set', columns='gene')
    step = step + 1
    cnt['step'] = step
    snp_qc = pd.concat([snp_qc, cnt])

    pos_imgt = vaf_normalize_to_normal(pos_imgt)

    ## Equalize number of SNPs for each allele
    pos_imgt = flip_snps(pos_imgt, vaf_col='obsVafNorm', new_col='obsVafNorm_flip')

    logging.info(f'SNP counts after filtering:\n{snp_qc}')
    ### Allele specific counts
    cnt = pos_imgt.allele.value_counts()
    cnt.name = 'snps'
    cnt.index.rename('allele', inplace=True)
    cnt = pd.DataFrame(cnt).sort_index()
    logging.info(f'\n{cnt}')
    
    ## Fit SNPs
    varsFile = f'{outDir}/{specimen_id}_workVars.tsv'
    pos_imgt.to_csv(varsFile, sep='\t', index=False)
    cmd = fitSequenza_cmd(
        specimen_id = specimen_id,
        varsFile = varsFile,
        outputDir = outDir,
        altMode = 'y',
        sequenzaModelRData = sequenzaModelRData,
        sequenzaTools = package_file_path(hlacn, 'sequenzaTools.R')
    )
    logging.info(' '.join(cmd))
    logfile = f'{logDir}/fitSequenza.log'
    handle = open(logfile, 'wt')
    job=sb.run(cmd, check=True, stdout=handle, stderr=sb.STDOUT)
    handle.close()

    # Copy the selected solution to the output directory for downstream analysis
    allelesFile = f'{outDir}/alt{args.solution}/{specimen_id}_alleles.tsv'
    if os.path.exists(allelesFile):
        logging.info(f'Copying {allelesFile} to {outDir} for downstream analysis')
        sb.run(['cp', '-f', allelesFile, outDir], check=True)
    else:
        raise FileNotFoundError(f'{allelesFile} does not exist')

    logging.info(f'{os.path.basename(__file__)} finished')
    logging.shutdown()

if __name__ == "__main__":
    main()
