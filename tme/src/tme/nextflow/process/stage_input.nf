/*
    Generate TME report
*/

include { saveFiles; getSoftwareName } from '../lib/functions'

params.options = [:]

process STAGE_INPUT {
 //   label 'process_medium'
    tag "${meta.id}"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: {
            filename -> saveFiles(filename:filename, options:[:], publish_dir:params.publish_dir)
            }

    input:
    tuple val(meta), path(manifest), path(input_folder)

    output:
    tuple val(meta), path("*.yml"), emit: tme_config

    script:
    """
    stage_input_tme.py ${manifest} ${input_folder}
    """
}
