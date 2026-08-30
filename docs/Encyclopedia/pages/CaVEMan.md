# CaVEMan

## Description ##

CaVEMan (Cancer Variants through Expectation Maximization) is a somatic variant calling algorithm used in the GDC for whole genome sequencing (WGS) data harmonization. It detects single nucleotide variants (SNVs) in tumor samples by comparing them to matched normal samples.

## Overview ##

CaVEMan is one of the four pipelines used for WGS variant calling at the GDC. Variant calling is performed with CaVEMan using tumor and normal alignments and generates single nucleotide variation (SNV) data.

### Input

* Tumor BAM - WGS
* Normal BAM - WGS

### Output

* VCF (Data Type: Raw Simple Somatic Mutation)

## References ##

1. [DNA-Seq Analysis Pipeline](/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/)
1. [Whole Genome Sequencing Variant Calling](/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/#whole-genome-sequencing-variant-calling)

## External Links ##

* [Overview of GDC Harmonization Workflows](https://github.com/NCI-GDC/gdc-workflow-overview/blob/master/README.md)
* [GDC Data Portal](https://portal.gdc.cancer.gov)

Categories: Workflow Type