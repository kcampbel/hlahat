/*
    Generate TME report
*/

include { saveFiles; getSoftwareName; generateTimestamp } from '../lib/functions'

params.options = [:]

process TME_REPORT {
 //   label 'process_medium'
    tag "${meta.id}"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: {
            filename -> saveFiles(filename:filename, options:[:], publish_dir:params.publish_dir)
            }
    conda params.conda_basedir + params.conda_envt

    input:
    tuple val(meta), path(config)
    tuple path(metadata), path(bctpm), path(counts)

    output:
    path "*.html"

    script:
    timestamp =  generateTimestamp()
    outfile = meta.study + '_' + meta.patient_id + '_' + meta.dob + '_BINF_' + meta.specimen_id + '_' + timestamp + '_tme_v0.1.html' 
    """
    cp -Lr "${params.rmd}" .
    echo '\
     rmarkdown::render( \
     input="R/main.Rmd", \
     output_file="../${outfile}", \
     params=list( \
      config_yml="../${config}", \
      metadata_f="../${metadata}", \
      counts_f="../${counts}", \
      bctpm_f="../${bctpm}") \
     )' > run.R

    echo '\
     rmarkdown::render( \
     input="R/main.Rmd", \
     output_file="../${outfile}", \
     params=list( \
      config_yml="../${config}", \
      metadata_f="../${metadata}", \
      counts_f="../${counts}", \
      bctpm_f="../${bctpm}") \
     )' >> run.R
    Rscript run.R
    """
}