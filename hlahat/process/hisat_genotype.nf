/*
    Extract HLA reads
*/

include { saveFiles; getSoftwareName } from '../modules/functions'

params.options = [:]

process EXTRACT_READS {
 //   label 'process_medium'
    tag "${meta}"
//    publishDir "${params.outdir}",
//        mode: params.publish_dir_mode,
//        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'merged_fastq', meta:meta, publish_by_meta:['id']) }
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: {
            filename -> saveFiles(filename:filename, options:[:], publish_dir:params.publish_dir)
            }

    input:
    tuple val(meta), path(reads1), path(reads2)
    path hisat_prefix

    output:
    tuple val(meta), path("*.extracted*.fq.gz") , emit: reads

    script:
    """
    export PATH=/opt/hisat2/hisat2-hisat2_v2.2.0_beta:/opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_scripts:/opt/samtools/bin:$PATH
    export PYTHONPATH=/opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_modules:$PYTHONPATH
    
    /opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_extract_reads_v_KC.py --base ${hisat_prefix}/genotype_genome \
        -p ${task.cpus} \
        -1 ${reads1} -2 ${reads2} \
        --database-list hla
    """
}

process HISAT_GENOTYPE {
 //   label 'process_medium'
    tag "${meta}"
//    publishDir "${params.outdir}",
//        mode: params.publish_dir_mode,
//        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'pipeline_info', meta:[:], publish_by_meta:[]) }
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: {
            filename -> saveFiles(filename:filename, options:[:], publish_dir:params.publish_dir)
            }

    input:
    tuple val(meta), path(reads)
    path hisat_prefix

    output:
    tuple val(meta), path("*.report") , emit: hla_types

    script:
    def software = getSoftwareName(task.process)
    """
    export PATH=/opt/hisat2/hisat2-hisat2_v2.2.0_beta:/opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_scripts:/opt/samtools/bin:$PATH
    export PYTHONPATH=/opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_modules:$PYTHONPATH
    
    /opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_locus_v_KC.py --genotype-genome ${hisat_prefix}/genotype_genome \
        -p ${task.cpus} \
        -1 ${reads[0]} -2 ${reads[1]} \
        --base hla --output-base ${meta} --keep-low-abundance-alleles > /dev/null
        #--base hla --locus-list A,B,C --output-base ${meta} --keep-low-abundance-alleles > /dev/null
    """
}