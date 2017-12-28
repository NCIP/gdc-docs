# Aligned Reads #
## Description ##

A read is a sequence obtained from a single sequencing experiment. An aligned read, is a sequence that has been aligned to a common reference genome. Typically these reads can number from the hundreds of thousands to tens of millions.

## Overview ##
The GDC supports the submission of aligned reads, in addition to unaligned reads. A data file containing aligned reads can be used as input for most GDC workflows<sup>1</sup>. During harmonization,  reads are aligned to the GRCh38 human genome with standardized protocols based on data type<sup>2,3</sup>. Generated aligned read files also contain unaligned reads to facilitate the retrieval of raw data by end users.

Aligned reads are available at the GDC Data Portal for:
* Whole Exome Sequencing
* Whole Genome Sequencing
* Transcriptome Sequencing

### Data Formats ###
Aligned reads are maintained in Binary Alignment Map (BAM) format.

## References ##
1. [GDC Data Dictionary - Submitted Aligned Reads](https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=submitted_aligned_reads)
2. [DNA-Seq Documentation](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/)
3. [RNA-Seq Documentation](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/)

## External Links ##
* N/A

Categories: Data Type
