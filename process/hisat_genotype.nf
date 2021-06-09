/*
    Extract HLA reads
*/

include { saveFiles; getSoftwareName } from '../modules/functions'

process EXTRACT_READS {
 //   label 'process_medium'
    tag "${meta.id}"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: {
            filename -> saveFiles(filename:filename, options:[:], publish_dir:params.publish_dir)
            }

    input:
    tuple val(meta), path(reads)
    path hisat_prefix

    output:
    //path "*.extracted*.fq.gz", emit: reads
    tuple val(meta), path("*.extracted*.fq.gz") , emit: reads

    script:
    def software = getSoftwareName(task.process)
    //def readList = reads.collect{ it.toString() }
    """
    export PATH=/opt/hisat2/hisat2-hisat2_v2.2.0_beta:/opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_scripts:/opt/samtools/bin:$PATH
    export PYTHONPATH=/opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_modules:$PYTHONPATH
    
    /opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_extract_reads_v_KC.py --base ${hisat_prefix}/genotype_genome \
        -1 ${reads[0]} -2 ${reads[1]} \
        --database-list hla
    """
}

process HISAT_GENOTYPE {
 //   label 'process_medium'
    tag "${meta.id}"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: {
            filename -> saveFiles(filename:filename, options:[:], publish_dir:params.publish_dir)
            }

    input:
    tuple val(meta), path(reads)
    path hisat_prefix
    //tuple val(), val(exp_id), path(R1), path(R2)

    output:
    tuple val(meta), path("*.report") , emit: hla_types

    script:
    def software = getSoftwareName(task.process)
    //def readList = reads.collect{ it.toString() }
    """
    export PATH=/opt/hisat2/hisat2-hisat2_v2.2.0_beta:/opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_scripts:/opt/samtools/bin:$PATH
    export PYTHONPATH=/opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_modules:$PYTHONPATH
    
    /opt/hisat2/hisat2-hisat2_v2.2.0_beta/hisatgenotype_locus_v_KC.py --genotype-genome ${hisat_prefix}/genotype_genome \
        -1 ${reads[0]} -2 ${reads[1]} \
        --base hla --locus-list A,B,C --output-base ${meta.id}.A,B,C --keep-low-abundance-alleles > /dev/null
    """
}