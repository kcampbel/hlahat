/*
    HLA typing compare tool
*/

include { saveFiles; getSoftwareName } from '../lib/functions'

params.options = [:]

process HLATYPE_COMPARE {
 //   label 'process_medium'
    tag "${meta.id}"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: {
            filename -> saveFiles(filename:filename, options:[:], publish_dir:params.publish_dir)
            }
    conda params.conda_basedir + params.conda_envt

    input:
    tuple val(meta), val(data), file(hisat2_hlatypes)
        
    output:
    tuple val(meta), val(data), file('*hlatypes_merged.tsv'), emit: hlatypes
    path ".command*"

    script:
    hlatypes_merged = meta.id + "_hlatypes_merged.tsv"
    """
    compare_hla_typing.py ${data.epic_hlatypes} ${hisat2_hlatypes} -o ${hlatypes_merged}
    """
}