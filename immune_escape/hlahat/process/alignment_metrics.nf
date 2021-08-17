/*
    Generate HLA-associated alignment metrics
*/
include { saveFiles; getSoftwareName } from '../lib/functions'

params.options = [:]

process ALIGNMENT_METRICS {
 //   label 'process_medium'
    tag "${meta.id}"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: {
            filename -> saveFiles(filename:filename, options:[:], publish_dir:params.publish_dir)
            }
    //conda params.conda_basedir + params.conda_envt

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta), path(reads), path(hla_fasta), path(hla_bed)

    output:
    tuple val(meta), path("*flagstat.txt"), path("*.fixpileups.txt"), emit: metrics
    path "*.pileups.txt"
    path "*.idxstats.txt"
    path ".command*"

    script:
    def software = getSoftwareName(task.process)
    """
    samtools flagstat ${bam} > ${meta.id}.flagstat.txt
    samtools idxstats ${bam} > ${meta.id}.idxstats.txt
    samtools mpileup ${bam} -l ${hla_bed} -q ${params.minMQ} -Q ${params.minBQ} > ${meta.id}.pileups.txt
    Rscript /code/summarize_pileups.R ${hla_bed} ${meta.id}.pileups.txt ${meta.id}.fixpileups
    """
}