/*
    Generate TME report
*/

include { saveFiles; getSoftwareName } from '../lib/functions'

params.options = [:]

process STAGE_INPUT {
 //   label 'process_medium'
    tag "${meta}"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: {
            filename -> saveFiles(filename:filename, options:[:], publish_dir:params.publish_dir)
            }
    conda params.conda_basedir + params.conda_envt

    input:
    tuple val(meta), path(manifest), path(input_folder)

    output:
    path '*.tsv', emit: tsv
    path '*.yml'
    path ".command*"

    script:
    """
    stage_input_tme.py ${manifest} ${input_folder}
    """
}
