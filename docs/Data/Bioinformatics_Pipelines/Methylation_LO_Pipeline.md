# Methylation Liftover Pipeline

## Introduction

The [Methylation Liftover Pipeline](https://gdc-docs.nci.nih.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=methylation_liftover_workflow) uses data from the Illumina Human Methylation 27 and 450 arrays to measure the level of methylation at known CpG sites. Likelihood and level of methylation are measured with beta values. Beta values are calculated from array intensities, which are associated with metadata such as chromosomal position, gene location, and distance from transcription start sites (TSS).

Beta values and probe IDs are lifted directly from existing TCGA level 3 methylation array data, which was processed using the hg19 reference genome. Existing coordinates for CpG sites were re-annotated to GRCh38 reference genome coordinates. Coordinates were then used to identify the associated transcripts, the associated CpG island (CGI), and the CpG sites' distance from these features.


## Methylation Beta Values Table Format


| Field | Definition |
|---|---|
| Composite Element | A unique ID for the array probe associated with a CpG site |
| Beta Value | Represents the ratio between the methylated array intensity and total array intensity, falls between 0 (lower levels of methylation) and 1 (higher levels of methylation) |
| Chromosome | The chromosome in which the probe binding site is located |
| Start | The start of the CpG site on the chromosome |
| End | The end of the CpG site on the chromosome |
| Gene Symbol | The symbol for genes associated with the CpG site  |
| Gene Type | A general classification of each gene (e.g. protein coding, miRNA, pseudogene) |
| Transcript ID |  Ensembl transcript IDs for each transcript associated with the CpG site |
| Position to TSS | Distance in base pairs from the CpG site to each associated transcript's start site  |
| CGI_Coordinate | The start and end coordinates of the CpG island associated with the CpG site |
| Feature_Type | The position of the CpG site in reference to the island: Island, N_Shore or S_Shore (0-2 kb upstream or downstream from CGI), or N_Shelf or S_Shelf (2-4 kb upstream or downstream from CGI) |


| I/O | Entity | Format |
|---|---|---|
| Input | [Submitted Methylation Beta Values](Data_Dictionary/viewer/#?view=table-definition-view&id=submitted_methylation_beta_value) |  TXT |
| Output | [Methylation Beta Values](Data_Dictionary/viewer/#?view=table-definition-view&id=methylation_beta_value) or Masked Copy Number Segment | TXT  |

[1] XXXXXXXXXX
