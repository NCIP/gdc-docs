# Variant Call Format (VCF) #

## Description ##
The Variant Call Format (VCF) is a standardized format for storing and reporting genomic sequence variations.[1]

## Overview ##
VCF files are used to report sequence variations (e.g., SNPs, indels and larger structural variants) together with rich annotations. VCF files are modular where the annotations and genotype information for a variant are separated from the call itself. VCF version 4.1 is the currently active format specification.

### Structure ###

VCF files are tab-delimited files that report a mutation for each row.  VCFs are available in a raw format or an annotated format that contains additional information about the location and consequences of each somatic variant. 

Details about the structure of the VCF is available in the [VCF 4.1 Specification](https://samtools.github.io/hts-specs/VCFv4.1.pdf). Changes made to the VCF format in support of the GDC are available in the [GDC VCF Format document](https://gdc-docs.nci.nih.gov/Data/File_Formats/VCF_Format/).

## References ##
1. [Variant Call Format](https://wiki.nci.nih.gov/display/TCGA/Variant+Call+Format)
2. [GDC VCF Format](https://gdc-docs.nci.nih.gov/Data/File_Formats/VCF_Format/)

## External Links ##
* [VCF 4.1 Specification](https://samtools.github.io/hts-specs/VCFv4.1.pdf)

Categories: Data Type
