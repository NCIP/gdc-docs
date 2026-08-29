# SeSAMe Methylation Beta Estimation

## Description ##

SeSAMe Methylation Beta Estimation is a workflow used in GDC methylation array harmonization.

## Overview ##

SeSAMe Methylation Beta Estimation is used for methylation array harmonization at the GDC. Methylation array harmonization is performed with SeSAMe using raw tumor or normal methylation array files and generates beta values and masked intensities data, which removes potential genotype information.

### Input

* Raw Methylation Array

### Output

* IDAT (Data Type: Masked Intensities)
* TXT (Data Type: Methylation Beta Value)

## References ##

1. [Methylation Analysis](/Data/Bioinformatics_Pipelines/Methylation_Pipeline/)
1. [SeSAMe Methylation Beta Values File Format](/Data/Bioinformatics_Pipelines/Methylation_Pipeline/#sesame-methylation-beta-values-file-format)

## External Links ##

* [SeSAMe - SEnsible Step-wise Analysis of Methylation data](https://github.com/zwdzwd/sesame)
* [GDC Data Portal](https://portal.gdc.cancer.gov)

Categories: Workflow Type