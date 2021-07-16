/*
    Mark Duplicates
*/
include { saveFiles; getSoftwareName } from '../lib/functions'

params.options = [:]

process PICARD_MARK_DUPLICATES {
 //   label 'process_medium'
    tag "${meta.id}"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: {
            filename -> saveFiles(filename:filename, options:[:], publish_dir:params.publish_dir)
            }

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.bam"), path("*.bam.bai"), emit:hla_bam_nodup

    script:
    def software = getSoftwareName(task.process)
    """
    picard MarkDuplicates I=${bam} O=${meta.id}.md.bam\
     ASSUME_SORT_ORDER=coordinate\
     METRICS_FILE=${meta.id}.md.txt\
     QUIET=true COMPRESSION_LEVEL=0\
     VALIDATION_STRINGENCY=LENIENT
    samtools index ${meta.id}.md.bam
    """
}