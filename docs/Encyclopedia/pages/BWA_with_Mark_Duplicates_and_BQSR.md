# BWA with Mark Duplicates and BQSR

## Description ##

BWA with Mark Duplicates and BQSR (Base Quality Score Recalibration) is a sequencing alignment and quality improvement workflow used in the GDC ATAC-Seq, targeted sequencing, whole genome sequencing (WGS) and whole exome sequencing (WXS) data harmonization. 

## Overview ##

BWA with Mark Duplicates and BQSR is a BWA based alignment workflow for ATAC-Seq, targeted sequencing, whole genome sequencing (WGS) and whole exome sequencing (WXS) with additional steps of Mark Duplicates and Base Quality Score Recalibration (BQSR). BWA with Mark Duplicates and BQSR is used for aligning sequence reads to a reference genome.

### Input

* Submitted sequencing reads

### Output

* BAM (Data Type: Aligned Reads)

## References ##

1. [DNA-Seq Analysis Pipeline](/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/)
1. [Alignment Workflow](/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/#alignment-workflow)

## External Links ##

* [Overview of GDC Harmonization Workflows](https://github.com/NCI-GDC/gdc-workflow-overview/blob/master/README.md)
* [BWA on GitHub](https://github.com/lh3/bwa)
* [GDC Data Portal](https://portal.gdc.cancer.gov)

Categories: Workflow Type