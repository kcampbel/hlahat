#!/usr/bin/env python
import sys
import os
import pandas as pd
import argparse
import subprocess as sb
import logging
from hlacn.lib.phasing import bcftools_cmd, bcf_to_df, fitSequenza_cmd

def parse_args(args=None):
    parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('specimen_id')
    parser.add_argument('normal_dna_bam')
    parser.add_argument('tumor_dna_bam')
    parser.add_argument('snp_bed')
    parser.add_argument('snp_tsv')
    parser.add_argument('seqz_file')
    parser.add_argument('sequenzaModelRData')
    parser.add_argument('outDir')
    
    return parser.parse_args(args)

def main():
    args = parse_args()
#    args = parse_args(
#    [
#        'PACT125_T_202824',
#        '/media/nfs/data/Workspace/Users/csmith/nf-GATK_Exome_Preprocess/asterand/PACT125_N_202825/PACT125_N_202825/PACT125_N_202825_normal_dna-final.bam',
#        '/media/nfs/data/Workspace/Users/csmith/nf-GATK_Exome_Preprocess/asterand/PACT125_T_202824/PACT125_T_202824/PACT125_T_202824_tumor_dna-final.bam',
#        '/home/csmith/csmith/hla/sequenza/snp2allele/vaf/snps/DNA_PACT125_N_202825_snps.bed',
#        '/home/csmith/csmith/hla/sequenza/snp2allele/vaf/snps/DNA_PACT125_N_202825_snps.tsv',
#        '/media/nfs/data/Workspace/Users/csmith/hla/sequenza/yma/PACT125_T_202824_gamma80_kmin13_thin/out/sequenza/PACT125_T_202824.hla.small.seqz.gz',
#        '/media/nfs/data/Workspace/Users/csmith/hla/sequenza/snp2allele/PACT125_T_202824_gamma80_kmin13_thin/out/sequenza/sequenzaModel.RData',
#        '/media/nfs/data/Workspace/Users/csmith/hla/sequenza/snp2allele/vaf/test',
#    ])
    specimen_id = args.specimen_id
    nBam = args.normal_dna_bam
    tBam = args.tumor_dna_bam
    posBed = args.snp_bed
    snp_f = args.snp_tsv
    seqzFile = args.seqz_file
    sequenzaModelRData = args.sequenzaModelRData
    outDir = f'{args.outDir}/{args.specimen_id}'
    logDir = f'{outDir}/log'

    missing = list()
    for ii in [tBam, nBam, posBed, snp_f, seqzFile, sequenzaModelRData]:
        if not os.path.exists(ii):
            missing.append(ii)
    if missing:
        raise FileNotFoundError(missing)

    genome_fa = '/media/nfs/data/References/hg19/hs37d5.chr/hs37d5.chr.fa'
    if not os.path.exists(logDir):
        os.makedirs(logDir, exist_ok=True)

    formatter = '%(asctime)s:%(levelname)s:%(name)s:%(funcName)s: %(message)s'
    logging.basicConfig(filename=f'{logDir}/snp2allele.log', format=formatter, filemode='w', level=logging.DEBUG)

    # BCFTOOLS
    posTsv = f'{outDir}/{specimen_id}_pos.tsv'
    cmd = bcftools_cmd(nBam, tBam, posBed, genome_fa, posTsv)
    logging.info(cmd)
    logfile = f'{logDir}/bcftools.log'
    handle = open(logfile, 'wt')
    job = sb.run(cmd, stdout=handle, stderr=sb.STDOUT, shell=True, check=True)
    handle.close()

    # Intersect with IMGT SNPs
    # Necessary step?
    pos_df = pd.read_table(posTsv, skiprows=1, names=['chromosome', 'position', 'ref', 'alt', 'nGT', 'tGT', 'nAD', 'tAD'])
    snp_df = pd.read_table(snp_f)
    snp_df = snp_df.rename(columns={'snp': 'alt', 'base.ref': 'ref'})
    snp_df['gene'] = snp_df.allele.str[0]

    tmp = bcf_to_df(pos_df, min_reads=0)
    pos_m = tmp.merge(snp_df, how='outer', indicator=True)
    logging.info(f'bcftools snps (left), IMGT snps (right)\n{pos_m._merge.value_counts()}')

    pos_imgt = pos_m[
    (pos_m._merge == 'both') &
    (pos_m.normalCounts > 20)
    ]

    # Drop shared SNPs
    snp_uniq = (pos_imgt.groupby(['position']).pos_imgt.count() == 1).reset_index()
    snp_uniq = snp_uniq[snp_uniq.pos_imgt]\
    .drop('pos_imgt', axis=1)
    pos_imgt = pos_imgt.merge(snp_uniq)
    pos_imgt = pos_imgt.astype({'start': 'int32', 'end': 'int32'})
    logging.info(f'Dropping SNPs shared between alleles\n{pos_imgt.gene.value_counts()}')

    # Intersect with Sequenza SNPs
    seqz = pd.read_table(seqzFile)
    seqz = seqz.rename(columns={'chr':'chromosome'})
    tmp = pos_imgt.drop(columns='_merge')
    imgt_seqz = tmp.merge(seqz[['chromosome', 'position']], how='outer', indicator=True)
    #imgt_seqz.groupby('gene')._merge.value_counts()
    imgt_seqz = imgt_seqz[imgt_seqz._merge == 'both']
    imgt_seqz = imgt_seqz.astype({'start': 'int32', 'end': 'int32'})
    logging.info(f'Intersect with Sequenza seqz.gz\n{imgt_seqz.gene.value_counts()}')

    # Fit SNPs
    varsFile = f'{outDir}/{specimen_id}_workVars.tsv'
    imgt_seqz.to_csv(varsFile, sep='\t', index=False)
    cmd = fitSequenza_cmd(
        specimen_id = specimen_id,
        varsFile = varsFile,
        outputDir = outDir,
        altMode = 'n',
        sequenzaModelRData = sequenzaModelRData,
        sequenzaTools = '/home/csmith/git/immune_escape/hlacn/lib/sequenzaTools.R'
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