# RNA-Seq #
## Description ##

RNA-Seq produces nucleotide sequence data as well as gene expression levels. RNA-Seq data originates from an RNA extract that was reverse transcribed and subsequently sequenced using a next-generation sequencing platform.  The number of reads determined to have originated from each transcript (usually
by alignment) are proportional to their expression level.

## Overview ##

The GDC harmonizes RNA-Seq data by aligning raw reads to the [GRCh38](LINK) reference genome build and calculating gene expression levels with standardized protocols. RNA-Seq data is mostly available for tumor samples, although RNA-Seq data is available for some normal samples.

### Data ###

RNA-Seq data is available as aligned reads (BAM), raw counts (TXT), and normalized (TXT) with [FPKM](LINK) or [FPKM-UQ](LINK). Reads that did not align are also included in BAM alignments to facilitate the retrieval of the original raw data.   

### Validation ###
### Analysis ###
## References ##
1. [GDC mRNA Expression Pipeline](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/)

## External Links ##
*

Categories: Experimental Strategy
