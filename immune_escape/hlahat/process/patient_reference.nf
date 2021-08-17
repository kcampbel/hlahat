/*
    Generate patient HLA reference
*/
include { saveFiles; getSoftwareName } from '../lib/functions'

params.options = [:]

process PATIENT_REFERENCE {
 //   label 'process_medium'
    tag "${meta.id}"
//    publishDir "${params.outdir}",
//        mode: params.publish_dir_mode,
//        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'pipeline_info', meta:[:], publish_by_meta:[]) }
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: {
            filename -> saveFiles(filename:filename, options:[:], publish_dir:params.publish_dir)
            }
    //conda params.conda_basedir + params.conda_envt

    input:
    tuple val(meta), path(hisatgt_report)
    path imgthla

    output:
    tuple val(meta.specimen_id), path("*custom_hla.fasta"), path("*custom_hla.allelic_differences.bed"), emit: patient_reference
    //tuple val(meta), path("*custom_hla.fasta"), path("*custom_hla.allelic_differences.bed"), emit: patient_reference
    tuple val(meta.specimen_id), path("*top_hlatypes.tsv") , emit: top_hlatypes
    //path "*hlatypes.tsv" 
    path ".command*"

    shell:
    def software = getSoftwareName(task.process)
    '''
    gen=`find !{imgthla}/msf/*gen.msf |tr '\\n' ','`
    nuc=`find !{imgthla}/msf/*nuc.msf |tr '\\n' ','`
    Rscript /code/generate_reference_files.R !{meta.specimen_id} !{hisatgt_report} !{params.nfields} ${gen} ${nuc}
    '''
}