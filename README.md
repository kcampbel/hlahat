# HLA-HAT: HLA Haplotype Analysis Toolkit

**Allele-specific analysis of the HLA genes from DNA and RNA sequencing data.**

📖 **Documentation:** https://hlahat.readthedocs.io

The human leukocyte antigen (HLA) genes are among the most polymorphic in the genome, which makes them notoriously hard to analyze against a linear, haploid reference: reads from divergent alleles mismap, and allele-specific signal is lost. **HLA-HAT** solves this by typing each sample's HLA alleles, building a **custom, personalized HLA reference**, and re-aligning reads to it — enabling accurate variant calling, allele-specific expression, and allelic-imbalance / loss-of-heterozygosity (LOH) analysis at the HLA locus.

This matters for cancer immunogenomics in particular: HLA LOH and allele-specific dysregulation are mechanisms of immune escape, and they are invisible to standard reference-based pipelines.

## How it works

HLA-HAT wraps [HISAT2](https://daehwankimlab.github.io/hisat2/) and [HISAT-genotype](http://ccb.jhu.edu/hisat-genotype/) and adds allele-specific reference construction and downstream analysis, in three modules:

1. **Preparation & assignment** — extract HLA-mapping reads and type each sample's HLA alleles.
2. **Reference generation** — report the HLA haplotype, build a personalized HLA reference FASTA, and realign HLA reads to it.
3. **Downstream analysis** — allele-specific alignment metrics, allelic-imbalance quantitation, and HLA variant detection, with dedicated workflows for **paired tumor–normal** and **tumor-only** designs.

## Installation & usage

HLA-HAT can be run four ways; see the [Usage guide](https://hlahat.readthedocs.io/en/latest/usage.html) for exact commands and parameters:

- **Docker** — the recommended path; the container bundles HISAT2, HISAT-genotype, and all dependencies.
- **WDL** — a workflow for scalable, reproducible runs.
- **Terra** — run the WDL on the cloud.
- **R** — the downstream analysis functions (metrics, allelic imbalance, visualization).

```bash
# Sketch — see the docs for the full command and required inputs
docker pull <hlahat-image>
docker run --rm -v "$PWD:/data" <hlahat-image> \
    hlahat --sample SAMPLE_ID --bam /data/sample.bam --out /data/out/
```

**Inputs:** aligned DNA/RNA sequencing reads (BAM), optionally paired tumor–normal.
**Outputs:** HLA types, a personalized HLA reference, allele-specific alignment metrics, allelic-imbalance estimates, and HLA variant calls.

## Documentation

Full walkthrough — introduction to HLA nomenclature, usage, and a module-by-module analysis guide — at **https://hlahat.readthedocs.io**.

## Citation & acknowledgment

Any academic or scholarly publication arising from the use of this software or derivative works should include the following acknowledgment:

> The software used in this research was created by Katie Campbell of UCLA.

# License
[GPL-3.0](https://github.com/kcampbel/hlahat/blob/main/LICENSE)

Any academic or scholarly publication arising from the use of this Software or any derivative works thereof will include the following acknowledgment: The Software used in this research was created by Katie Campbell of UCLA.

