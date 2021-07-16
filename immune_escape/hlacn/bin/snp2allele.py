#!/usr/bin/env python
import sys
import os
import pandas as pd
import argparse
import subprocess as sb
import logging
from hlacn.lib.command import bcftools_cmd, fitSequenza_cmd
from hlacn.lib.munge import bcf_to_df, vaf_normalize_to_normal, flip_snps

def parse_args(args=None):
    parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--specimen_id', required=True, help='Specimen id')
    parser.add_argument('-n', '--normal_dna_bam', required=True, help='Normal BAM')
    parser.add_argument('-t', '--tumor_dna_bam', required=True, help='Tumor BAM')
    parser.add_argument('-g', '--genome_fasta', required=True, help='Reference genome FASTA')
    parser.add_argument('-s', '--snp_tsv', required=True, help='IMGT SNPs tsv')
    parser.add_argument('-m', '--sequenzaModelRData', required=True, help='Sequenza model RData file')
    parser.add_argument('-o', '--outDir', required=True, help='Output dir')
    
    return parser.parse_args(args)

def main():
    args = parse_args()
#    args = parse_args(
#    [
#        '-i', 'PACT125_T_202824',
#        '-n', '/media/nfs/data/Workspace/Users/csmith/nf-GATK_Exome_Preprocess/asterand/PACT125_N_202825/PACT125_N_202825/PACT125_N_202825_normal_dna-final.bam',
#        '-t' ,'/media/nfs/data/Workspace/Users/csmith/nf-GATK_Exome_Preprocess/asterand/PACT125_T_202824/PACT125_T_202824/PACT125_T_202824_tumor_dna-final.bam',
#        '-g', '/media/nfs/data/References/hg19/hs37d5.chr/hs37d5.chr.fa',
#        '-s', '/home/csmith/csmith/hla/sequenza/snp2allele/vaf/snps/DNA_PACT125_N_202825_snps.tsv',
#        '-m', '/media/nfs/data/Workspace/Users/csmith/hla/sequenza/snp2allele/PACT125_T_202824_gamma80_kmin13_thin/out/sequenza/sequenzaModel.RData',
#        '-o', '/media/nfs/data/Workspace/Users/csmith/hla/sequenza/snp2allele/vaf/test/snp_count_ratio_flipped',
#    ])
    specimen_id = args.specimen_id
    nBam = args.normal_dna_bam
    tBam = args.tumor_dna_bam
    snp_f = args.snp_tsv
    genome_fasta = args.genome_fasta
    sequenzaModelRData = args.sequenzaModelRData
    outDir = f'{args.outDir}/{args.specimen_id}'
    logDir = f'{outDir}/log'

    missing = list()
    for ii in [tBam, nBam, snp_f, sequenzaModelRData]:
        if not os.path.exists(ii):
            missing.append(ii)
    if missing:
        raise FileNotFoundError(missing)

    if not os.path.exists(logDir):
        os.makedirs(logDir, exist_ok=True)

    formatter = '%(asctime)s:%(levelname)s:%(name)s:%(funcName)s: %(message)s'
    logging.basicConfig(filename=f'{logDir}/snp2allele.log', format=formatter, filemode='w', level=logging.DEBUG)

    # BCFTOOLS
    ## Make SNP bed
    snp_df = pd.read_table(snp_f)
    snp_df = snp_df.rename(columns={'snp': 'alt', 'base.ref': 'ref'})
    snp_bed = snp_df.assign(start=snp_df.position-1).loc[:, ['chromosome', 'start', 'position']]
    snp_bed_f = f'{outDir}/{specimen_id}_snps.bed'
    snp_bed.to_csv(snp_bed_f, sep='\t', index=False, header=False)

    posTsv = f'{outDir}/{specimen_id}_pos.tsv'
    cmd = bcftools_cmd(nBam, tBam, snp_bed_f, genome_fasta, posTsv)
    logging.info(cmd)
    logfile = f'{logDir}/bcftools.log'
    handle = open(logfile, 'wt')
    job = sb.run(cmd, stdout=handle, stderr=sb.STDOUT, shell=True, check=True)
    handle.close()

    # Intersect with IMGT SNPs
    pos_df = pd.read_table(posTsv, skiprows=1, names=['chromosome', 'position', 'ref', 'alt', 'nGT', 'tGT', 'nAD', 'tAD'])

    tmp = bcf_to_df(pos_df, min_reads=0)
    pos_m = tmp.merge(snp_df, how='outer', indicator=True)
    logging.info(f'bcftools snps (left), IMGT snps (right)\n{pos_m._merge.value_counts()}')

    # Filter SNPs
    min_normal = 20
    pos_imgt = pos_m[pos_m._merge == 'both']
    nsnps = pos_imgt.shape[0]
    pos_imgt = pos_imgt[pos_imgt.normalCounts > min_normal]
    logging.info(f'Removing {nsnps - pos_imgt.shape[0]} SNPs with less than {min_normal} reads in normal' + 
    f'\n{pos_imgt.gene.value_counts()}')
    
    # Drop shared SNPs between allele 1 and allele 2
    snp_uniq = (pos_imgt.groupby(['position']).pos_imgt.count() == 1).reset_index()
    snp_uniq = snp_uniq[snp_uniq.pos_imgt].drop('pos_imgt', axis=1)
    pos_imgt = pos_imgt.merge(snp_uniq)
    pos_imgt = pos_imgt.astype({'start': 'int32', 'end': 'int32'})
    logging.info(f'Dropping SNPs shared between alleles\n{pos_imgt.gene.value_counts()}')

    # Normalize tumor VAF to normal
    pos_imgt = vaf_normalize_to_normal(pos_imgt)

    # Equalize number of SNPs for each allele
    pos_imgt = flip_snps(pos_imgt, vaf_col='obsVafNorm', new_col='obsVafNorm_flip')

    # Fit SNPs
    varsFile = f'{outDir}/{specimen_id}_workVars.tsv'
    pos_imgt.to_csv(varsFile, sep='\t', index=False)
    cmd = fitSequenza_cmd(
        specimen_id = specimen_id,
        varsFile = varsFile,
        outputDir = outDir,
        altMode = 'n',
        sequenzaModelRData = sequenzaModelRData,
        sequenzaTools = f'{os.path.dirname(__file__)}/../lib/sequenzaTools.R'
    )
    logging.info(' '.join(cmd))
    logfile = f'{logDir}/fitSequenza.log'
    handle = open(logfile, 'wt')
    job=sb.run(cmd, check=True, stdout=handle, stderr=sb.STDOUT)
    handle.close()
    allelesFile = f'{outDir}/alt1/{specimen_id}_alleles.tsv'
    sb.run(['ln', '-fs', allelesFile, outDir], check=True)

    logging.info(f'{specimen_id} finished.')
    logging.shutdown()

if __name__ == "__main__":
    main()
