# PACT RNAseq expression matrix update and batch correction pipeline

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
  * genes.tsv sample RNAseq transcriptome quantification file

```bash
conda activate process_exprs

process_expression.py /path/to/manifest /path/to/epic/output folder
```

## Outputs

The pipeline will update:
* PACT expression matrix with hgnc symbols
* PACT expression matrix with gene ids 
* PACT expression matrix logfile
* PACT/Xena batch corrected TPM expression matrix
* PACT/Xena batch corrected TPM logfile
* PACT/Xena metadata file

The location of input and output paths for the matrices, logs, and metadata are configured in `src/process_exprs/data/pipeline.yml`. A custom `pipeline.yaml` can be specified with `process_expression.py`.

