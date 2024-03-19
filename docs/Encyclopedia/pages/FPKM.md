# FPKM #
## Description ##

Fragments Per Kilobase of transcript per Million mapped reads (FPKM) is a simple expression level normalization method. The FPKM normalizes read count based on gene length and the total number of mapped reads.  

## Overview ##

FPKM is implemented at the GDC on gene-level read counts that are produced by STAR<sup>1</sup> and generated using custom scripts<sup>2</sup>. The formula used to generate FPKM values is as follows:

FPKM = [RM<sub>g</sub> * 10<sup>9</sup> ] / [RM<sub>t</sub> * L]

* RM<sub>g</sub>: The number of reads mapped to the gene
* RM<sub>t</sub>: The total number of read mapped to protein-coding sequences in the alignment
* L: The length of the gene in base pairs

The scalar (10<sup>9</sup>) is added to normalize the data to "__kilo__ base" and "__million__ mapped reads."

FPKM files are available as tab delimited files with the Ensembl gene IDs in the first column and the expression values in the second. See [FPKM-UQ](FPKM-UQ.md) for an alternative method of gene expression level normalization.


## References ##
1. [STAR-Fusion pipeline](/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/#star-fusion-pipeline)
2. [GDC mRNA-Seq Documentation](/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/)


## External Links ##
* [Ensembl Human Genome](http://www.ensembl.org/Homo_sapiens/Info/Annotation)
* [GENCODE 36](https://www.gencodegenes.org/human/release_36.html)

Categories: Workflow Type
