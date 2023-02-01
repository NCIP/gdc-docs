# RNA-Seq #
## Description ##

RNA-Seq is a sequencing method used to determine gene expression levels. RNA-Seq data originates from extracted RNA that was reverse transcribed into DNA and sequenced on a next-generation sequencing platform. The number of reads determined to have originated from each transcript (usually by alignment) are proportional to their expression level.

## Overview ##

RNA-Seq data in the GDC is used to generate a gene expression profile for tumor samples across many cancer types and to determine which gene expression levels are responsible for tumor development. The GDC harmonizes RNA-Seq data by aligning raw RNA reads to the GRCh38 reference genome build and calculating gene expression levels with standardized protocols<sup>1</sup>. RNA-Seq data is mostly available for tumor samples, although some normal samples have associated RNA-Seq data.

### Data ###

RNA-Seq data is available as aligned reads (BAM) and expression levels as: raw counts and normalized with TPM, FPKM, or FPKM-UQ. Reads that did not align are also included in BAM files to facilitate the retrieval of the original raw data.   

## References ##
1. [GDC mRNA Expression Pipeline](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/)

## External Links ##
* N/A

Categories: Experimental Strategy
