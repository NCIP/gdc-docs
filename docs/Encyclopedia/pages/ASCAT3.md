# ASCAT3

## Description ##

ASCAT3 is an advanced copy number variation (CNV) analysis pipeline used in GDC genotyping array harmonization.

## Overview ##

ASCAT3 improves upon previous versions to provide more accurate CNV detection in tumor and normal samples from genotyping arrays. It processes the raw array files to generate high-resolution segmented CNV data.

### Input

* Tumor CEL - Genotyping Array
* Normal CEL - Genotyping Array

### Output

* TXT (Data Type: Allele-specific Copy Number Segment) 
* TSV (Data Type: Gene Level Copy Number)

## References ##

1. [Copy Number Variation Analysis Pipeline](/Data/Bioinformatics_Pipelines/CNV_Pipeline/)
1. [ASCAT Pipelines](/Data/Bioinformatics_Pipelines/CNV_Pipeline/#ascat-pipelines)

## External Links ##

* [Vanloo lab](https://github.com/VanLoo-lab/ascat/tree/master/ReleasedData/TCGA_SNP6_hg38)
* [GDC Data Portal](https://portal.gdc.cancer.gov)

Categories: Workflow Type