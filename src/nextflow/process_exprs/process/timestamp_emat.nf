/*
    Update PACT expression set, log and metadata
*/

include { saveFiles; getSoftwareName } from '../lib/functions'

params.options = [:]

process TIMESTAMP_EMAT {
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
    tuple val(meta), path(exprs), path(exprs_log), path(metadata_tsv)
    val(outdir)

    output:
    path ".command*"

    script:
    """
    update_timestamp.py\
     --exprs ${exprs}\
     --exprs_log ${exprs_log}\
     --metadata_tsv ${metadata_tsv}\
     -o ${outdir}
    """
}
