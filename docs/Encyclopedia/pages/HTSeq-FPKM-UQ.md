# HTSeq-FPKM-UQ #
## Description ##

Fragments Per Kilobase of transcript per Million mapped reads upper quartile (FPKM-UQ) is a RNA-Seq-based expression normalization method.  The FPKM-UQ is based on a modified version of the [FPKM](HTSeq-FPKM.md) normalization method.  

## Overview ##

FPKM-UQ is implemented at the GDC on gene-level read counts that are produced by HTSeq<sup>1</sup> and generated using custom scripts<sup>2</sup>. The formula used to generate FPKM-UQ values is as follows:

FPKM = [RM<sub>g</sub> * 10<sup>9</sup> ] / [RM<sub>75</sub> * L]

* RM<sub>g</sub>: The number of reads mapped to the gene
* RM<sub>75</sub>: The number of read mapped to the 75th percentile gene in the alignment.
* L: The length of the gene in base pairs

Like [HTSeq - count](HTSeq-Counts.md) files, FPKM-UQ files are available as tab delimited files with the Ensembl gene IDs in the first column and the expression values in the second.

### Notes
- The scalar (10<sup>9</sup>) is added to normalize the values to "__kilo__ base" and "__million__ mapped reads."
- FPKM-UQ values tend to be much higher than FPKM values because of the large difference between the total mapped number of reads in an alignment and the mapped number of reads to one gene.  

### Tools ###
1. [HTSeq Website](http://htseq.readthedocs.io/)

## References ##
1. Anders, S., Pyl, P.T. and Huber, W., 2014. HTSeq-a Python framework to work with high-throughput sequencing data. Bioinformatics, p.btu638.
2. [GDC mRNA-Seq Documentation](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/)


## External Links##
* [Ensembl Human Genome](http://www.ensembl.org/Homo_sapiens/Info/Annotation)
* [GENCODE 22](https://www.gencodegenes.org/human/release_22.html)

Categories: Workflow Type
