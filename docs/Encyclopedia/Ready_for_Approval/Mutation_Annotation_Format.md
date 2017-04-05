# Mutation Annotation Format (MAF) #
## Description ##
A Mutation Annotation Format (MAF) is a tab-delimited file containing somatic and/or germline mutation annotations. MAF files containing any germline mutation annotations are kept in the controlled access portion of the GDC Data Portal. MAF files containing only somatic mutations are kept in the open access portion of the GDC Data Portal.[1]

## Overview ##
Variants are discovered by aligning DNA sequences derived from tumor samples and sequences derived from normal samples to a reference sequence. A MAF file identifies, for each sample, the discovered putative or validated mutations and categorizes those mutations (polymorphism, deletion, or insertion) as somatic (originating in the tumor tissue) or germline (originating from the germline) as well as the annotation for those mutations.[1]

MAF files are generated at the GDC using the Variant Aggregation pipeline.

### Structure ###
The structure of the MAF is available in the [GDC MAF Specification] (https://gdc-docs.nci.nih.gov/Data/File_Formats/MAF_Format/) along with descriptions for each field.

## References ##
1. [GDC MAF Format](https://gdc-docs.nci.nih.gov/Data/File_Formats/MAF_Format/)
2. [Mutation Annotation Format](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format)

## External Links ##
* [TCGA MAF specification](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification)


Categories: Data Format
