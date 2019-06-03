# Methylation Liftover Pipeline

## Introduction

The [DNA Methylation Liftover Pipeline](https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=methylation_liftover_workflow) uses data from the Illumina Infinium Human Methylation 27 (HM27) and HumanMethylation450 (HM450) arrays to measure the level of methylation at known CpG sites as beta values, calculated from array intensities (Level 2 data) as Beta = M/(M+U).

Using probe sequence information provided in the manufacturer's manifest, HM27 and HM450 probes were remapped to the GRCh38 reference genome [1]. Type II probes with a mapping quality of <10, or Type I probes for which the methylated and unmethylated probes map to different locations in the genome, and/or had a mapping quality of <10, had an entry of '\*' for the 'chr' field, and '-1' for coordinates. These coordinates were then used to identify the associated transcripts from GENCODE v22, the associated CpG island (CGI), and the CpG sites' distance from each of these features. Multiple transcripts overlapping the target CpG were separated with semicolons. Beta values were inherited from existing TCGA Level 3 DNA methylation data (hg19-based) based on Probe IDs.


## Methylation Beta Values Table Format

Descriptions for fields present in GDC Harmonized Methylation Beta Values Table are detailed below:

| Field | Definition |
|---|---|
| Composite Element | A unique ID for the array probe associated with a CpG site |
| Beta Value | Represents the ratio between the methylated array intensity and total array intensity, falls between 0 (lower levels of methylation) and 1 (higher levels of methylation) |
| Chromosome | The chromosome in which the probe binding site is located |
| Start | The start of the CpG site on the chromosome |
| End | The end of the CpG site on the chromosome |
| Gene Symbol | The symbol for genes associated with the CpG site. Genes that fall within 1,500 bp upstream of the transcription start site (TSS) to the end of the gene body are used.    |
| Gene Type | A general classification for each gene (e.g. protein coding, miRNA, pseudogene) |
| Transcript ID |  Ensembl transcript IDs for each transcript associated with the genes detailed above |
| Position to TSS | Distance in base pairs from the CpG site to each associated transcript's start site  |
| CGI Coordinate | The start and end coordinates of the CpG island associated with the CpG site |
| Feature Type | The position of the CpG site in reference to the island: Island, N_Shore or S_Shore (0-2 kb upstream or downstream from CGI), or N_Shelf or S_Shelf (2-4 kbp upstream or downstream from CGI) |

---
| I/O | Entity | Format |
|---|---|---|
| Input | [Submitted Methylation Beta Values](/Data_Dictionary/viewer/#?view=table-definition-view&id=submitted_methylation_beta_value) |  TXT |
| Output | [Methylation Beta Values](/Data_Dictionary/viewer/#?view=table-definition-view&id=methylation_beta_value) or Masked Copy Number Segment | TXT  |

## File Access and Availability

| Type | Description | Format |
|---|---|---|
| Methylation Beta Value | A table that associates array probes with CpG sites and associated metadata. |  TXT |


[1]. Zhou, Wanding, Laird Peter L., and Hui Shen. "Comprehensive characterization, annotation and innovative use of Infinium DNA methylation BeadChip probes." Nucleic Acids Research. (2016): doi: 10.1093/nar/gkw967
