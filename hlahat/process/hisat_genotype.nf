/*
    HISAT GENOTYPE
*/

include { saveFiles; getSoftwareName } from '../modules/functions'

params.options = [:]

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
    tuple val(meta), path("*.report") , emit: hisatgt_report

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