# HTSeq-Counts #
## Description ##

HTSeq is a Python package that calculates the number of mapped reads to each protein-coding gene.

## Overview ##

The first step in generating gene expression values from an RNA-Seq alignment at the GDC is generating a count of the reads mapped to each gene. These counts are performed using HTSeq and are calculated at the gene level. HTSeq-Count files are available in a tab-delimited format with one Ensembl gene ID column and one mapped reads column for each gene. These files are then processed further with custom scripts to generate [FPKM](HTSeq-FPKM.md) and [FPKM-UQ](HTSeq-FPKM-UQ.md) values.

### Tools ###
1. [HTSeq Website](http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html)

## References ##
1. [GDC mRNA-Seq Documentation](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/)
2. Anders, S., Pyl, P.T. and Huber, W., 2014. HTSeq–a Python framework to work with high-throughput sequencing data. Bioinformatics, p.btu638.

Categories: Workflow Type
