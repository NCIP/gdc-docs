# SNP Array-Based Data #
## Description ##

SNP Arrays produce data based on DNA hybridization levels to a set of array probes. Hybridization occurs on the arrays based on the presence of probe-specific SNPs in the chromosome.

## Overview ##

SNP Array-Based data is used at the GDC to perform Copy Number Variation analyses. A copy number segmentation analysis is performed using DNACopy. Chromosomal segments are identified and their copy number are estimated based on similarities between adjacent SNP array probe binding sites<sup>1</sup>.   

### Data Formats ###

The GDC stores processed SNP Array-Based data in the active harmonized portal. They are stored as tab-delimited copy number segmentation files, which identify the segment regions and their estimated copy number. The raw array intensity files in CEL format along with downstream copy number variation and genotyping analyses are available in the GDC Legacy Archive<sup>2</sup>.  

## References ##
1. [GDC Copy Number Variation Documentation](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/CNV_Pipeline/)
2. [GDC Legacy Archive](https://portal.gdc.cancer.gov/legacy-archive/)

## External Links ##
* TBD

Categories: Data Type
