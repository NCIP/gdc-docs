# SNP Array-Based Data #
## Description ##

SNP Arrays produce data based on DNA hybridization levels to a set of array probes. Hybridization occurs on the arrays based on the presence or absence of probe-specific SNPs in the chromosome.

## Overview ##

SNP Array-Based data is used at the GDC to perform Copy Number Variation analyses. A copy number segmentation analysis is performed using DNACopy. Chromosomal segments are identified and their copy numbers are estimated based on similarities between adjacent SNP array probe binding sites<sup>1</sup>.   

### Data Formats ###

The GDC stores processed SNP Array-Based data in the active harmonized portal. They are stored as tab-delimited copy number segmentation files, which identify the segment regions and their estimated copy number. The raw array intensity files in CEL format are available in the GDC Data Portal, and the downstream copy number variation and genotyping analyses are available as Birdseed files.

## References ##
1. [GDC Copy Number Variation Documentation](/Data/Bioinformatics_Pipelines/CNV_Pipeline/)

## External Links ##
* N/A

Categories: Data Type
