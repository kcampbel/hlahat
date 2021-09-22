.. hlahat documentation master file, created by
   sphinx-quickstart on Fri May 28 12:49:24 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

HLA-HAT: HLA Haplotype Analysis Toolkit
==================================
The HLA Haplotype Analysis Toolkit (**HLA-HAT**) was designed to support analysis of the human leukocyte antigen (HLA) genes in the allele-specific context. The polymorphic nature of these genes makes it especially difficult to analyze and visualize allele-specific information in the context of the linear, haploid reference genome. Thus, HLA-HAT implements `HISAT2`_ and `HISAT-genotype`_ for the extraction and typing from HLA reads and constructs a custom reference HLA genome for downstream variant calling, allelic imbalance quantitation, and additional analysis.

The User Guide
==================================

This documentation describes how to get started with HLA-HAT, including the extraction, assignment, and downstream analysis of HLA allele-specific sequencing reads.

.. toctree::
   :maxdepth: 3

   intro
   usage
   analysis
   resources
   references

.. _`HISAT2`: https://daehwankimlab.github.io/hisat2/
.. _`HISAT-genotype`: http://ccb.jhu.edu/hisat-genotype/index.php/Main_Page
