# Tumor Microenvironment Reporting Pipeline

## Installation

Dependencies:
* conda

```
bash install.sh
```

## Running the pipeline 

### Local run

EPIC inputs:
* manifest.yml
* EPIC outputs folder including
  * segments.txt from cellularity-with-sequenza
  * Annotated.tsv from epitope-annotation

External inputs:
* Batch corrected TPM matrix
* PACT/Xena metadata tsv
* PACT raw counts matrix 

```
conda activate tme

sid='PACT280C_T_PP001128'
manifest=inputs/PACT280C/PACT280C_T_PP001128/manifest.yml
input_dir=inputs/PACT280C/PACT280C_T_PP001128

bctpm='/media/nfs/data/References/xena/subsets/PactTcgaTargetGtex_rsem_gene_tpm/exprs_combat_Breast_hgnc.tsv.gz'
counts="${HOME}/git/bioinfo-fio/tme/test_data/pact_rsem_raw.exprs.gid.tsv.gz"
metadata="${HOME}/git/bioinfo-fio/tme/test_data/meta_pact_xena.tsv"

tme_report.py ${manifest} ${input_dir} -e "'--bctpm,${bctpm},--counts,${counts},--metadata,${metadata}'" 
```

## Outputs

The pipeline will output a html report in the following format:

```
{study id}_{patient id}_{patient dob}_BINF_{specimen id}_{timestamp}_tme_{tme report version}.html`
```

