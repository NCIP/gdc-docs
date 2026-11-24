# BWA-aln #
## Description ##

BWA-aln is an alignment algorithm that is part of the Burrows-Wheeler Alignment Tool used to align short reads to a large reference genome.

## Overview ##

BWA-aln is implemented in the GDC DNA-Seq and miRNA-Seq pipelines and used for the initial alignment of reads to the reference genome. BWA-aln is used in the DNA-Seq pipeline only if the mean read length is less than 70 bp, whereas [BWA-mem](LINK) is used if the mean read length is greater than or equal to 70 bp.    

### Tools ###
## References ##
1. Li, H. and Durbin, R., 2010. Fast and accurate long-read alignment with Burrows-Wheeler transform. Bioinformatics, 26(5), pp.589-595.

## External Links ##
* [GDC DNA-Seq Documentation](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/)

Categories: Workflow Type
