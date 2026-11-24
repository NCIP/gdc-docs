# STAR - Counts

## Description ##

STAR - Counts is an RNA expression pipeline used in GDC RNA-Seq harmonization.

## Overview ##

STAR - Counts is a pipeline used for to quantify gene expression from RNA-Seq data. Quantification is performed with STAR using tumor or normal alignments and generates transcriptome profiling data.

### Input

* Genomic BAM - RNA-Seq

### Output

* TSV (Data Type: Gene Expression Quantification, Splice Junction Quantification)

[![STAR Pipeline Diagram](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/images/RNA-Seq-DR32_Image.png)](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/images/RNA-Seq-DR32_Image.png "Click to see the full image.")

## References ##

1. [mRNA Analysis Pipeline](/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/)
1. [mRNA Expression Workflow](/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/#mrna-expression-workflow)
1. [mRNA Quantification Command Line Parameters](/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/#mrna-quantification-command-line-parameters)

## External Links ##

* [STAR](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
* [Overview of GDC Harmonization Workflows](https://github.com/NCI-GDC/gdc-workflow-overview/blob/master/README.md)
* [GDC Data Portal](https://portal.gdc.cancer.gov)

Categories: Workflow Type