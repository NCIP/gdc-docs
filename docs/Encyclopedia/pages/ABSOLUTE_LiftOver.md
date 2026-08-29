# ABSOLUTE LiftOver

## Description ##

ABSOLUTE LiftOver is a copy number variation (CNV) pipeline used in GDC genotyping array harmonization.

## Overview ##

The ABSOLUTE LiftOver workflow in the GDC is a copy number variation (CNV) pipeline used for genotyping array harmonization. It converts genomic coordinates from hg19 to GRCh38, ensuring high-quality, curated CNV data for downstream analysis. This pipeline also provides gene-level copy numbers, along with purity and ploidy measurements.


### Input

* Tumor CEL - Genotyping Array

### Output

* TSV (Data Type: Gene Level Copy Number)

## References ##

1. [Copy Number Variation Analysis Pipeline](/Data/Bioinformatics_Pipelines/CNV_Pipeline/)
1. [ABSOLUTE Copy Number](/Data/Bioinformatics_Pipelines/CNV_Pipeline/#absolute-copy-number)

## External Links ##

* [TCGA PanCancer analysis papers](https://doi.org/10.1016/j.ccell.2018.03.007)
* [GDC Data Portal](https://portal.gdc.cancer.gov)

Categories: Workflow Type