/*
    Generate patient HLA reference
*/

include { saveFiles; getSoftwareName } from '../lib/functions'

params.options = [:]

process ALIGN_HLA_READS {
 //   label 'process_medium'
    tag "${meta.id}"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: {
            filename -> saveFiles(filename:filename, options:[:], publish_dir:params.publish_dir)
            }

    input:
    tuple val(meta), path(reads), path(hla_fasta), path(hla_bed)
    //tuple val(meta), path(reads)
    //tuple path(hla_fasta), path(hla_bed)

    output:
    tuple val(meta), path("*.bam"), path("*.bam.bai"), emit:hla_bam

    script:
    def software = getSoftwareName(task.process)
    if ( meta.seqtype.equals('tumor_rna') )
    """
    /opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisat2-build ${hla_fasta} custom_hla.hisat
    /opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisat2 \
        --no-discordant \
        -p ${task.cpus} \
        -1 ${reads[0]} -2 ${reads[1]} \
        -S ${meta.id}.sam \
        -x custom_hla.hisat 
    /opt/samtools/bin/samtools sort -@ ${task.cpus} -o ${meta.id}.bam ${meta.id}.sam
    /opt/samtools/bin/samtools index ${meta.id}.bam
    """
    else
    """
    /opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisat2-build ${hla_fasta} custom_hla.hisat
    /opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisat2 \
        --no-spliced-alignment \
        --no-discordant \
        -p ${task.cpus} \
        -1 ${reads[0]} -2 ${reads[1]} \
        -S ${meta.id}.sam \
        -x custom_hla.hisat
    /opt/samtools/bin/samtools sort -@ ${task.cpus} -o ${meta.id}.bam ${meta.id}.sam
    /opt/samtools/bin/samtools index ${meta.id}.bam
    """
}