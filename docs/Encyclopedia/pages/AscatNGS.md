# AscatNGS

## Description ##

ASCATNGS is a copy number variation (CNV) analysis pipeline tailored for next-generation sequencing (NGS) data used in GDC whole genome sequencing (WGS) harmonization.

## Overview ##

ASCATNGS (__now deprecated at the GDC__) is specifically designed to detect CNVs in tumor and normal samples from NGS data. It processes aligned BAM files to generate segmented CNV data, providing high-resolution insights into genomic alterations.

### Input

* Tumor BAM - WGS
* Normal BAM - WGS

### Output

* TSV (Data Type: Gene Level Copy Number)
* TXT (Data Type: Copy Number Segment)

## References ##

1. [Copy Number Variation Analysis Pipeline](/Data/Bioinformatics_Pipelines/CNV_Pipeline/)
1. [ASCAT Pipelines](/Data/Bioinformatics_Pipelines/CNV_Pipeline/#ascat-pipelines)
1. [Whole Genome Sequencing Variant Calling](/Data/Bioinformatics_Pipelines/DNA_Seq_WGS/)

## External Links ##

* [AscatNGS](https://github.com/cancerit/ascatNgs)
* [GDC Data Portal](https://portal.gdc.cancer.gov)

Categories: Workflow Type