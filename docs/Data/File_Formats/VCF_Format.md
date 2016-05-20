# GDC VCF Format

## Introduction

The GDC VCF file format follows standards of the [Variant Call Format (VCF) Version 4.1 Specification](https://samtools.github.io/hts-specs/VCFv4.1.pdf). Raw Simple Somatic Mutation VCF files are unannotated, whereas Annotated Somatic Mutation VCF files include extensive, consistent, and pipeline-agnostic annotation of somatic variants.

## VCF file structure

A standard VCF file is composed of three parts, in the following order:  

### Part 1: Meta-information

The meta-information section of a VCF file comprises lines that begin with `##`. Some key components of this section are:

*  **INDIVIDUAL:** information on the study participant, including:
 	* *NAME:* Submitter ID (barcode) associated with the participant, and
	* *ID:* GDC case UUID
*  **SAMPLE:** sample information, including:
	*  *ID:* NORMAL or TUMOR
	*  *NAME:* Submitter ID (barcode) of the aliquot
	*  *ALIQUOT_ID:* GDC aliquot UUID
	*  *BAM_ID:* BAM file UUID
*   **FILTER:** description of filters that have been applied to the data
*   **FORMAT:** description of genotype fields
*   **INFO:** format of *additional information* fields (see [GDC INFO Fields](#gdc-info-fields) below)

### Part 2: Column Header Line

This is a single tab-delimited line that names the data columns in the following order:

0. **#CHROM:** chromosome
0. **POS:** position
0. **ID:** identifier
0. **REF:** reference base(s)
0. **ALT:** alternate base(s)
0. **QUAL:** quality
0. **FILTER:** filter status
0. **INFO:** additional information
0. **FORMAT:** format of sample genotype data
0. **NORMAL:** normal sample genotype data
0. **TUMOR:** tumor sample genotype data

See [Variant Call Format (VCF) Version 4.1 Specification](https://samtools.github.io/hts-specs/VCFv4.1.pdf) for details.

### Part 3: Data

This part of a VCF file contains tab-delimited information, one line per record, with each record representing a position in the genome.

## GDC INFO fields

The following variant annotation fields are currently included in Annotated Somatic Mutation VCF files.  Please refer to [GDC Somatic Mutation Annotation Pipeline](../Bioinformatics_Pipelines/Annotation_Pipeline.md) for details on how this information is generated.

|field|description|
|--|--|
|Allele|the variant allele used to calculate the consequence|
|Consequence|consequence type of this variant|
|IMPACT|the impact modifier for the consequence type|
|SYMBOL|the gene symbol|
|Gene|Ensembl stable ID of affected gene|
|Feature_type|type of feature. Currently one of Transcript, RegulatoryFeature, MotifFeature.|
|Feature|Ensembl stable ID of feature|
|BIOTYPE|Biotype of transcript or regulatory feature|
|EXON|the exon number (out of total number)|
|INTRON|the intron number (out of total number)|
|HGVSc|the HGVS coding sequence name|
|HGVSp|the HGVS protein sequence name|
|cDNA_position|relative position of base pair in cDNA sequence|
|CDS_position|relative position of base pair in coding sequence|
|Protein_position|relative position of amino acid in protein|
|Amino_acids|change in amino acids (only given if the variant affects the protein-coding sequence)|
|Codon|the alternative codons with the variant base in upper case|
|Existing_variation|known identifier of existing variant|
|ALLELE_NUM|Allele number from input; 0 is reference, 1 is first alternate etc|
|DISTANCE|Shortest distance from variant to transcript|
|STRAND|the DNA strand (1 or -1) on which the transcript/feature lies|
|FLAGS|transcript quality flags|
|VARIANT_CLASS|Sequence Ontology variant class|
|SYMBOL_SOURCE|the source of the gene symbol|
|HGNC_ID|HGNC gene ID|
|CANONICAL|a flag indicating if the transcript is denoted as the canonical transcript for this gene|
|TSL|Transcript support level|
APPRIS  
|CCDS|the CCDS identifer for this transcript, where applicable|
|ENSP|the Ensembl protein identifier of the affected transcript|
|SWISSPROT|UniProtKB/Swiss-Prot identifier of protein product|
|TREMBL|UniProtKB/TrEMBL identifier of protein product|
|UNIPARC|UniParc identifier of protein product|
|RefSeq|RefSeq gene ID|
|GENE_PHENO|Indicates if overlapped gene is associated with a phenotype, disease or trait|
|SIFT|the SIFT prediction and/or score, with both given as prediction(score)|
|PolyPhen|the PolyPhen prediction and/or score|
|DOMAINS|the source and identifer of any overlapping protein domains|
|HGVS_OFFSET|Indicates by how many bases the HGVS notations for this variant have been shifted|
|GMAF|Non-reference allele and frequency of existing variant in 1000 Genomes|
|AFR_MAF|Non-reference allele and frequency of existing variant in 1000 Genomes combined African population|
|AMR_MAF|Non-reference allele and frequency of existing variant in 1000 Genomes combined American population|
|ASN_MAF|Non-reference allele and frequency of existing variant in 1000 Genomes combined Asian population|
|EUR_MAF|Non-reference allele and frequency of existing variant in 1000 Genomes combined European population|
|EAS_MAF|Non-reference allele and frequency of existing variant in 1000 Genomes combined East Asian population|
|SAS_MAF|Non-reference allele and frequency of existing variant in 1000 Genomes combined South Asian population|
|AA_MAF|Non-reference allele and frequency of existing variant in NHLBI-ESP African American population|
|EA_MAF|Non-reference allele and frequency of existing variant in NHLBI-ESP European American population|
|ExAC_MAF|Frequency of existing variant in ExAC combined population|
|ExAC_Adj_MAF|Adjusted frequency of existing variant in ExAC combined population|
|ExAC_AFR_MAF|Frequency of existing variant in ExAC African/American population|
|ExAC_AMR_MAF|Frequency of existing variant in ExAC American population|
|ExAC_EAS_MAF|Frequency of existing variant in ExAC East Asian population|
|ExAC_FIN_MAF|Frequency of existing variant in ExAC Finnish population|
|ExAC_NFE_MAF|Frequency of existing variant in ExAC Non-Finnish European population|
|ExAC_OTH_MAF|Frequency of existing variant in ExAC combined other combined populations|
|ExAC_SAS_MAF|Frequency of existing variant in ExAC South Asian population|
|CLIN_SIG|Clinical significance of variant from dbSNP|
|SOMATIC|Somatic status of existing variant(s)|
|PHENO|Indicates if existing variant is associated with a phenotype, disease or trait|
|MOTIF_NAME|the source and identifier of a transcription factor binding profile aligned at this position|
|MOTIF_POS|The relative position of the variation in the aligned TFBP|
|HIGH_INF_POS|a flag indicating if the variant falls in a high information position of a transcription factor binding profile (TFBP)|
|MOTIF_SCORE_CHANGE|The difference in motif score of the reference and variant sequences for the TFBP|
|ENTREZ|Entrez ID|
|EVIDENCE|evidence of variant exist|
