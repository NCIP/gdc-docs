# STAR 2-Pass #
## Description ##

The STAR aligner is an alignment algorithm that specializes in RNA-Seq.  The STAR 2-Pass method first aligns the reads and identifies splice junctions and then uses the splice-junction information in a second pass to generate a higher-quality alignment.  

## Overview ##

RNA-Seq reads harmonized by the GDC are initially aligned to the reference genome using the STAR 2-Pass method. The resulting [BAM files](LINK) are used for downstream [gene expression quantification](LINK).

### Tools ###
## References ##
1. [GDC mRNA-Seq Documentation](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/)
2. [STAR github](https://github.com/alexdobin/STAR)

## External Links ##
* TBD

Categories: Workflow Type
