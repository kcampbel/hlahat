/*
    Extract HLA reads
*/
include { saveFiles; getSoftwareName } from '../lib/functions'

params.options = [:]

process EXTRACT_READS {
 //   label 'process_medium'
    tag "${meta.id}"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: {
            filename -> saveFiles(filename:filename, options:[:], publish_dir:params.publish_dir)
            }
    //conda params.conda_basedir + params.conda_envt

    input:
    tuple val(meta), path(reads1), path(reads2)
    path hisat_prefix

    output:
    tuple val(meta), path("*.extracted*.fq.gz"), emit: reads
    path ".command*"

    script:
    def software = getSoftwareName(task.process)
    if ( meta.seqtype.equals('tumor_rna') )
    """
    /opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_extract_reads_v_KC.py --base ${hisat_prefix}/genotype_genome \
        -p ${task.cpus} \
        -1 ${reads1} -2 ${reads2} \
        --is-rna \
        --database-list hla
    """
    else
    """
    /opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_extract_reads_v_KC.py --base ${hisat_prefix}/genotype_genome \
        -p ${task.cpus} \
        -1 ${reads1} -2 ${reads2} \
        --database-list hla
    """
}