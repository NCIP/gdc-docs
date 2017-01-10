# miRNA Expression Quantification #
## Description ##
A table that associates miRNA IDs with read count and a normalized count in reads-per-million-miRNA-mapped.
## Overview ##
Following alignment, BAM files are processed through the [miRNA Expression Workflow](https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=mirna_expression_workflow "miRNA Expression Workflow").
The outputs of the miRNA profiling pipeline report raw read counts and counts normalized to reads per million mapped reads (RPM) in two separate files mirnas.quantification.txt and isoforms.quantification.txt. The former contains summed expression for all reads aligned to known miRNAs in the miRBase reference. If there are multiple alignments to different miRNAs or different regions of the same miRNA, the read is flagged as cross-mapped and every miRNA annotation is preserved.
### Data Formats ###
TXT file
## References ##
1.https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/miRNA_Pipeline/

## External Links ##
* 

Categories: Data Type