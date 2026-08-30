# GATK4 MuTect2 Tumor-Only

## Description ##

GATK4 MuTect2 Tumor-Only is a somatic variant calling pipeline utilized in the GDC for analyzing tumor-only whole exome sequencing (WXS) and targeted sequencing data without normal samples. 

## Overview ##

GATK4 MuTect2 Tumor-Only is a pipeline used for WXS and targeted sequencing somatic variant calling at the GDC. Somatic variant calling is performed with GATK4 MuTect2 Tumor-Only using tumor alignments and generates single nucleotide variants (SNVs) data.

### Input

* Tumor BAM - WXS or Targeted Sequencing

### Output

* VCF (Data Type: Raw Simple Somatic Mutation)

## References ##

1. [DNA-Seq Analysis Pipeline](/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/)
1. [Tumor-Only Variant Calling Workflow](/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/#tumor-only-variant-calling-workflow)

## External Links ##

* [Overview of GDC Harmonization Workflows](https://github.com/NCI-GDC/gdc-workflow-overview/blob/master/README.md)
* [GDC Data Portal](https://portal.gdc.cancer.gov)

Categories: Workflow Type