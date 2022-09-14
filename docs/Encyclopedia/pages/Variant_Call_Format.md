# Variant Call Format (VCF) #

## Description ##
The Variant Call Format (VCF) is a standardized format for storing and reporting genomic sequence variations.

## Overview ##
VCF files are used to report sequence variations (e.g., SNPs, indels and larger structural variants) together with rich annotations. VCF files are modular where the annotations and genotype information for a variant are separated from the call itself. VCF version 4.1 is the currently active format specification<sup>1</sup>.

VCF files are generated at the GDC with one of four variant callers (MuSe, MuTect2, Pindel, VarScan) by comparing a tumor alignment to a normal alignment from the same patient. All GDC VCFs from patients are protected data and require dbGaP credentials to access.

VCF files are produced on a case-level for each variant caller (four VCF files per case). All VCF files within a project that were produced by a single pipeline are aggregated to produce a MAF file.  

### Structure ###

VCF files are tab-delimited files that report a mutation for each row. VCFs are available at the GDC in a raw format or an annotated format that contains additional information about the location and consequences of each somatic variant.

Details about the structure of the VCF is available in the VCF 4.1 Specification<sup>1</sup>. Changes made to the VCF format in support of the GDC are available in the GDC VCF Format documentation<sup>2</sup>.

## References ##
1. [VCF 4.1 Specification](https://samtools.github.io/hts-specs/VCFv4.1.pdf)
2. [GDC VCF Format](https://docs.gdc.cancer.gov/Data/File_Formats/VCF_Format/)

## External Links ##
* [VCF 4.1 Specification](https://samtools.github.io/hts-specs/VCFv4.1.pdf)

Categories: Data Type
