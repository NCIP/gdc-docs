# VarScan2

## Description ##

VarScan2 is a somatic variant calling pipeline used in GDC whole exome sequencing (WXS) and targeted sequencing harmonization.

## Overview ##

VarScan2 is one of the four pipelines used for WXS and targeted sequencing somatic variant calling at the GDC. Somatic variant calling is performed with VarScan2 using tumor and normal alignments and generates single-nucleotide polymorphism (SNP) data.

### Input

* Tumor BAM - WXS or Targeted Sequencing
* Normal BAM - WXS or Targeted Sequencing

### Output

* VCF (Data Type: Raw Simple Somatic Mutation)

[![VarScan2 Pipeline Diagram](https://gdc.cancer.gov/system/files/public/image/varscan-somatic-variant-calling-pipeline.png)](https://gdc.cancer.gov/system/files/public/image/varscan-somatic-variant-calling-pipeline.png "Click to see the full image.")

## References ##

1. [DNA Seq Processing at the GDC](/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/)
1. [Somatic Variant Calling Workflow](/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/#somatic-variant-calling-workflow)

## External Links ##

* [VarScan](https://dkoboldt.github.io/varscan/)
* [VarScan 2: somatic mutation and copy number alteration discovery in cancer by exome sequencing.](http://genome.cshlp.org/content/22/3/568.short)
* [Overview of GDC Harmonization Workflows](https://github.com/NCI-GDC/gdc-workflow-overview/blob/master/README.md)
* [GDC Data Portal](https://portal.gdc.cancer.gov)

Categories: Workflow Type