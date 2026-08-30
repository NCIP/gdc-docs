# MuSE

## Description ##

MuSE is a somatic variant calling pipeline used in GDC whole exome sequencing (WXS) and targeted sequencing harmonization.

## Overview ##

MuSE is one of the four pipelines used for WXS and targeted sequencing somatic variant calling at the GDC. Somatic variant calling is performed with MuSE using tumor and normal alignments and generates single-nucleotide polymorphism (SNP) data.

### Input
* Tumor BAM - WXS or Targeted Sequencing
* Normal BAM - WXS or Targeted Sequencing

### Output

* VCF (Data Type: Raw Simple Somatic Mutation)

[![MuSE Pipeline Diagram](https://gdc.cancer.gov/files/public/image/muse-somatic-variant-calling-pipeline.png)](https://gdc.cancer.gov/files/public/image/muse-somatic-variant-calling-pipeline.png "Click to see the full image.")

## References ##

1. [MuSE Command Line Parameters at the GDC](/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/#muse)
1. [DNA Seq Processing at the GDC](/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/)

## External Links ##

* [MuSE at MDACC](https://bioinformatics.mdanderson.org/public-software/muse/)
* [Overview of GDC Harmonization Workflows](https://github.com/NCI-GDC/gdc-workflow-overview/blob/master/README.md)
* [GDC Data Portal](https://portal.gdc.cancer.gov)

Categories: Workflow Type