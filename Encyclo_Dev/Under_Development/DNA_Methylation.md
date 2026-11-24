# DNA Methylation #

## Description ##
DNA methylation is an epigenetic mark which can be associated with transcriptional inactivity when located in promoter regions. In the context of the GDC, DNA methylation refers to the covalent modification of cytosine bases at the C-5 position, generally within a CpG sequence context.
## Overview ##
DNA methylation data contains information on raw and normalized signal intensities, detection confidence and calculated beta values for methylated and unmethylated probes. Higher level data includes masking of data, including SNP-related data. Details on data masking can be found in the DESCRIPTION.txt files that accompany data level and MAGE-TAB archives.
### Data ###
DNA methylation data consists of multiple files associated with each platform.

| Platform Code | Platform Alias | Platform Name | Platform File Description | Platform File Type | Platform File Example |
| --- | --- | --- | --- | --- | --- | --- |
| HumanMethylation27 | HumanMethylation27 | Illumina Infinium Human DNA Methylation 27 | Intensity data file with statistics for each bead type in terms of bead count, mean and standard deviation per dye | Binary (.idat) | TBD |
| HumanMethylation27 | HumanMethylation27 | Illumina Infinium Human DNA Methylation 27 | Calculated beta values and mean signal intensities for replicate methylated and unmethylated probes | Tab-delimited, ASCII text (.txt) | TBD |
| HumanMethylation27 | HumanMethylation27 | Illumina Infinium Human DNA Methylation 27 | Calculated beta values, gene symbols, chromosomes and genomic coordinates (build 36). Some data have been masked (including known SNPs)| Tab-delimited, ASCII text (.txt) | TBD |
| HumanMethylation450 | HumanMethylation450 | Illumina Infinium Human DNA Methylation 450 | Intensity data file with statistics for each bead type in terms of bead count, mean and standard deviation per dye | Binary (.idat) | TBD |
| HumanMethylation450 | HumanMethylation450 | Illumina Infinium Human DNA Methylation 450 | Background-corrected methylated (M) and unmethylated (U) summary intensities as extracted by the methylumi package | Tab-delimited, ASCII text (.txt) | TBD |
| HumanMethylation450 | HumanMethylation450 | Illumina Infinium Human DNA Methylation 450 | Calculated beta values, gene symbols, chromosomes and genomic coordinates (hg18). Some data have been masked (including known SNPs) | Tab-delimited, ASCII text (.txt) | |
| IlluminaDNAMethylation_OMA002_CPI | IlluminaDNAMethylation | Illumina DNA Methylation OMA002 Cancer Panel I | Cy3 and Cy5 signals and detection confidence of methylated probes | Tab-delimited, ASCII text (.txt) | TBD |
| IlluminaDNAMethylation_OMA002_CPI | IlluminaDNAMethylation | Illumina DNA Methylation OMA002 Cancer Panel I | Calculated beta values | Tab-delimited, ASCII text (.txt) | TBD |
| IlluminaDNAMethylation_OMA003_CPI | IlluminaDNAMethylation | Illumina DNA Methylation OMA003 Cancer Panel I | Cy3 and Cy5 signals and detection confidence of methylated probes | Tab-delimited, ASCII text (.txt) | TBD |
| IlluminaDNAMethylation_OMA003_CPI | IlluminaDNAMethylation | Illumina DNA Methylation OMA003 Cancer Panel I | Calculated beta values | Tab-delimited, ASCII text (.txt) | TBD |
### Validation ###
DNA Methylation data undergoes standard validation including:
- MD5 Validation

### Analysis ###
Analyzing DNA methylation data involves comparing ..

## References ##
1. [DNA Methylation](https://wiki.nci.nih.gov/display/TCGA/DNA+methylation "DNA Methylation")

## External Links ##
* N/A

Categories: Data Category
