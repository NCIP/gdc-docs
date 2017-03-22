# GDC MAF Format
- Version 1.0.0

## Introduction

Mutation Annotation Format (MAF) is a tab-delimited text file that describes mutations and are generated on a project-level.  The GDC produces MAF files at two permission levels: __protected__ and __somatic__ (or open-access). One MAF files is produced per variant calling pipeline per GDC project. MAFs are produced by aggregating the GDC annotated VCF files.

Annotated VCF files often have variants reported on multiple transcripts whereas the protected MAF file generated from the VCFs (\*protected.maf) only reports the most critically affected one. Somatic MAFs (\*somatic.maf) are further processed to remove low quality and potential germline variants. For tumor samples that contain variants from multiple combinations of tumor normal aliquot pairs, only one pair is selected in the Somatic MAF. Somatic MAFs are publicly available and can be freely distributed.

The GDC MAF file format is based on the [TCGA Mutation Annotation Format](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification) specifications, but some additional columns were included and some were removed.



### Important updates from TCGA v2.4

1. Column #4 "NCBI_Build" is GRCh38 by default
2. Column #60 VEP name "STRAND" is changed to "TRANSCRIPT_STRAND" to avoid confusion with Column#8 "Strand"
3. Column #15 "dbSNP_Val_Status" is not implemented
4. Column #26 "Mutation_Status": multiallelic calls from MuSE are labeled in this column as  MuSEMulti
5. Column #32 "Sequencer" includes sequencers used. If different sequencers were used to generate normal and tumor data, normal sequencer is listed first.
6. Column #113-116 "vcf_info", "vcf_format", "vcf_tumor_gt" and "vcf_normal_gt" are the corresponding columns from VCF files. Including them facilitates parsing specific variant information.
7. Column #109 "FILTER" entries are separated by both commas and semi-colons.
8. Column #117 "GDC_Validation_Status": GDC also collects TCGA validation sequences.  It compares these with variants derived from Next-Generation Sequencing data from the same sample and populates the comparison result in "GDC_Validation_Status".
    * "Valid", if Alternative allele(s) in the tumor validation sequence is(are) the same as GDC variant call
    * "Invalid", if none of the alternative allele(s) in the tumor validation sequence is the same as GDC variant call
    * "Inconclusive" if two alternative allele exists, and one matches while the other does not.
    * "Unknown" if no validation sequence exists
9. Column #118 GDC_Valid_Somatic is TRUE if GDC_Validation_Status is "Valid" and the variant is "Somatic" in validation calls.  It is FALSE if these criteria are not met.



## Somatic MAF file structure

This is very similar to the protected MAF format, but with some important additional filters.  The process for modifying a protected MAF into a somatic MAF are the following:

1. Aliquot Selection: only one tumor-normal pair are selected for each tumor sample based on the plate number, sample_type, analyte_type and other features extracted from tumor TCGA aliquot barcode.
2. Low quality variants filtering: currently variants with FILTER="PASS" are kept and n_depth < 8 are removed.
3. Germline Masking
    * Remove variants with Mutation_Status does not equal "Somatic"
    * Keep variants if any of the follow three criteria applies, but dbSNP sites that are not annotated as "Somatic" in VEP are also removed
        * GDC_Valid_Somatic is TRUE
        * VEP "Somatic" field is not empty
        * The variant affects transcript sequences. This criterion will be met if Variant_Classification is not 3'UTR, 3'Flank, 5'UTR, 5'Flank, IGR, or Intron and Variant_Classification is Frame_Shift_Del, Frame_Shift_Ins, In_Frame_Del, In_Frame_Ins, Missense_Mutation, Nonsense_Mutation, Silent, Splice_Site, Translation_Start_Site, Nonstop_Mutation, RNA, Targeted_Region, De_novo_Start_InFrame, or De_novo_Start_OutOfFrame.
4. Removal of the following columns:
    * vcf_info
    * vcf_format
    * vcf_tumor_gt
    * vcf_normal_gt
5. Set values to be blank in the following columns that may contain information about germline genotypes:
    * Match_Norm_Seq_Allele1
    * Match_Norm_Validation_Allele1
    * Match_Norm_Validation_Allele2
    * n_ref_count
    * n_alt_count

Note: the criteria for allowing mutations into open-access are purposefully implemented to overcompensate and filter out any germline mutations. If omission of true-positive somatic mutations is a concern, the GDC recommends using protected MAFs.  

## Protected MAF File Structure

The table below describes the columns in a protected MAF and their definitions. Note that the somatic (open-access) MAF structure is the same except for having the last four columns removed.

| Column | Description |
|---|---|
| 1 - Hugo_Symbol | [HUGO](http://www.genenames.org/) symbol for the gene (HUGO symbols are always in all caps). "Unknown" is used for regions that do not correspond to a gene. |
| 2 - Entrez_Gene_Id |Entrez gene ID (an integer).  "0" is used for regions that do not correspond to a gene|
| 3 - Center|One or more genome sequencing center reporting the variant.|
| 4 - NCBI_Build|Any TGCA accepted genome identifier. (e.g., hg18, hg19, GRCh37, GRCh38, GRCh37-lite, 36, 36.1, 37, 38)|
| 5 - Chromosome | The affected chromosome number |
| 6 - Start_Position | Lowest numeric position of the reported variant on the genomic reference sequence. Mutation start coordinate.|
| 7 - End_Position | Highest numeric genomic position of the reported variant on the genomic reference sequence. Mutation end coordinate.|
| 8 - Strand | Genomic strand of the reported allele. Currently, all variants will report the positive strand: '+'.|
| 9 - Variant_Classification | Translational effect of variant allele. |
| 10 - Variant_Type |Type of mutation. TNP (tri-nucleotide polymorphism) is analogous to DNP (di-nucleotide polymorphism) but for three consecutive nucleotides. ONP (oligo-nucleotide polymorphism) is analogous to TNP but for consecutive runs of four or more. (SNP, DNP, TNP, ONP, INS, DEL, or Consolidated)|
| 11 - Reference_Allele | The plus strand reference allele at this position. Include the deleted sequence for a deletion, or "-" for an insertion. |
| 12 - Tumor_Seq_Allele1 | Primary data genotype for tumor sequencing (discovery) allele 1. A "-" symbol for a deletion represents a variant. A "-" symbol for an insertion represents wild-type allele. Novel inserted sequence for insertion does not include flanking reference bases. |
| 13 - Tumor_Seq_Allele2 | Tumor sequencing (discovery) allele 2. |
| 14 - dbSNP_RS|The rs IDs from the dbSNP database, null if not found in any database used, or "novel" if there is no dbSNP record, but it is found in other databases|
| 15 - dbSNP_Val_Status| The dbSNP validation status. Semicolon-separated list of validation statuses. The union of all rs IDs is taken when there are multiple.|
| 16 - Tumor_Sample_Barcode| Aliquot barcode for the tumor sample, which includes the two additional fields indicating plate and well position.|
| 17 - Matched_Norm_Sample_Barcode| Aliquot barcode for the matched normal sample including the two additional fields indicating plate and well position.|
| 18 - Match_Norm_Seq_Allele1 | Primary data genotype. Matched normal sequencing allele 1.  A "-" symbol for a deletion represents a variant. A "-" symbol for an insertion represents wild-type allele. Novel inserted sequence for insertion does not include flanking reference bases. (cleared in somatic MAF) |
| 19 - Match_Norm_Seq_Allele2| Matched normal sequencing allele 2 |
| 20 - Tumor_Validation_Allele1| Secondary data from orthogonal technology. Tumor genotyping (validation) for allele 1. A "-" symbol for a deletion represents a variant. A "-" symbol for an insertion represents wild-type allele. Novel inserted sequence for insertion does not include flanking reference bases. |
| 21 - Tumor_Validation_Allele2| Secondary data from orthogonal technology. Tumor genotyping (validation) for allele 2.|
| 22 - Match_Norm_Validation_Allele1 | Secondary data from orthogonal technology. Matched normal genotyping (validation) for allele 1. A "-" symbol for a deletion represents a variant. A "-" symbol for an insertion represents wild-type allele. Novel inserted sequence for insertion does not include flanking reference bases (cleared in somatic MAF)|
| 23 - Match_Norm_Validation_Allele2| Secondary data from orthogonal technology. Matched normal genotyping (validation) for allele 2 (cleared in somatic MAF)|
| 24 - Verification_Status|Second pass results from independent attempt using same methods as primary data source. Generally reserved for 3730 Sanger Sequencing.|
| 25 - Validation_Status|Second pass results from orthogonal technology.|
| 26 - Mutation_Status|Updated to reflect validation or verification status and to be in agreement with the VCF VLS field. The values allowed in this field are constrained by the value in the Validation_Status field.|
| 27 - Sequencing_Phase|TCGA sequencing phase. Phase should change under any circumstance that the targets under consideration change.|
| 28 - Sequence_Source|Molecular assay type used to produce the analytes used for sequencing. Allowed values are a subset of the SRA 1.5 library_strategy field values. This subset matches those used at CGHub.|
| 29 - Validation_Method | The assay platforms used for the validation call|
| 30 - Score | Not in use.|
| 31 - BAM_File | Not in use.|
| 32 - Sequencer | Instrument used to produce primary sequence data|
| 33 - Tumor_Sample_UUID | BCR aliquot UUID for tumor sample|
| 34 - Matched_Norm_Sample_UUID | BCR aliquot UUID for matched normal sample|
| 35 - HGVSc | The coding sequence of the variant in HGVS recommended format|
| 36 - HGVSp | The protein sequence of the variant in HGVS recommended format|
| 37 - HGVSp_Short | Same as the HGVSp column, but using 1-letter amino-acid codes|
| 38 - Transcript_ID|The transcript that the variant affects|
| 39 - Exon_Number | The exon number (out of total number)|
| 40 - t_depth | Read depth across this locus in tumor BAM|
| 41 - t_ref_count | Read depth supporting the reference allele in tumor BAM|
| 42 - t_alt_count | Read depth supporting the variant allele in tumor BAM|
| 43 - n_depth | Read depth across this locus in normal BAM|
| 44 - n_ref_count | Read depth supporting the reference allele in normal BAM (cleared in somatic MAF)|
| 45 - n_alt_count | Read depth supporting the variant allele in normal BAM (cleared in somatic MAF)|
| 46 - all_effects | A semicolon delimited list of all possible variant effects, sorted by priority ([SYMBOL,Consequence,HGVSp_Short,Transcript_ID,RefSeq])|
|47 - Allele|The variant allele used to calculate the consequence|
|48 - Gene|Stable Ensembl ID of affected gene|
|49 - Feature|Stable Ensembl ID of feature (transcript, regulatory, motif)|
|50 - Feature_type|type of feature. Currently one of Transcript, RegulatoryFeature, MotifFeature (or blank)|
|51 - Consequence | Consequence type of this variation; sequence ontology terms |
|52 - cDNA_position | Relative position of base pair in the cDNA sequence |
|53 - CDS_position | Relative position of base pair in coding sequence |
|54 - Protein_position | Relative position of affected amino acid in protein |
|55 - Amino_acids | Only given if the variation affects the protein-coding sequence |
|56 - Codons | The alternative codons with the variant base in upper case|
|57 - Existing_variation | Known identifier of existing variation|
|58 - ALLELE_NUM | Allele number from input; 0 is reference, 1 is first alternate etc|
|59 - DISTANCE | Shortest distance from variant to transcript|
|60 - TRANSCRIPT_STRAND | The DNA strand (1 or -1) on which the transcript/feature lies|
|61 - SYMBOL | The gene symbol|
|62 - SYMBOL_SOURCE | The source of the gene symbol|
|63 - HGNC_ID | Gene identifier from the HUGO Gene Nomenclature Committee|
|64 - BIOTYPE | Biotype of transcript|
|65 - CANONICAL | A flag (YES) indicating that the VEP-based canonical transcript, the longest translation, was used for this gene |
|66 - CCDS | The CCDS identifier for this transcript, where applicable|
|67 - ENSP | The Ensembl protein identifier of the affected transcript|
|68 - SWISSPROT | UniProtKB/Swiss-Prot accession|
|69 - TREMBL | UniProtKB/TrEMBL identifier of protein product|
|70 - UNIPARC | UniParc identifier of protein product|
|71 - RefSeq | RefSeq identifier for this transcript|
|72 - SIFT | The SIFT prediction and/or score, with both given as prediction (score)|
|73 - PolyPhen | The PolyPhen prediction and/or score|
|74 - EXON | The exon number (out of total number)|
|75 - INTRON | The intron number (out of total number)|
|76 - DOMAINS | The source and identifier of any overlapping protein domains|
|77 - GMAF | Non-reference allele and frequency of existing variant in [1000 Genomes](http://www.internationalgenome.org/)|
|78 - AFR_MAF | Non-reference allele and frequency of existing variant in 1000 Genomes combined African population|
|79 - AMR_MAF | Non-reference allele and frequency of existing variant in 1000 Genomes combined American population|
|80 - ASN_MAF | Non-reference allele and frequency of existing variant in 1000 Genomes combined Asian population|
|81 - EAS_MAF | Non-reference allele and frequency of existing variant in 1000 Genomes combined East Asian population|
|82 - EUR_MAF | Non-reference allele and frequency of existing variant in 1000 Genomes combined European population|
|83 - SAS_MAF | Non-reference allele and frequency of existing variant in 1000 Genomes combined South Asian population|
|84 - AA_MAF | Non-reference allele and frequency of existing variant in [NHLBI-ESP](http://evs.gs.washington.edu/EVS/) African American population|
|85 - EA_MAF | Non-reference allele and frequency of existing variant in NHLBI-ESP European American population|
|86 - CLIN_SIG | Clinical significance of variant from dbSNP|
|87 - SOMATIC |Somatic status of each ID reported under Existing_variation (0, 1, or null)|
|88 - PUBMED | Pubmed ID(s) of publications that cite existing variant|
|89 - MOTIF_NAME |The source and identifier of a transcription factor binding profile aligned at this position|
|90 - MOTIF_POS | The relative position of the variation in the aligned TFBP|
|91 - HIGH_INF_POS | A flag indicating if the variant falls in a high information position of a transcription factor binding profile (TFBP) (Y, N, or null)|
| 92 - MOTIF_SCORE_CHANGE | The difference in motif score of the reference and variant sequences for the TFBP |
| 93 - IMPACT | The impact modifier for the consequence type|
| 94 - PICK | Indicates if this block of consequence data was picked by VEP's pick feature (1 or null)|
| 95 - VARIANT_CLASS | Sequence Ontology variant class |
| 96 - TSL | Transcript support level |
| 97 - HGVS_OFFSET | Indicates by how many bases the HGVS notations for this variant have been shifted |
| 98 - PHENO|Indicates if existing variant is associated with a phenotype, disease or trait (0, 1, or null) |
| 99 - MINIMISED | Alleles in this variant have been converted to minimal representation before consequence calculation (1 or null) |
| 100 - ExAC_AF | Global Allele Frequency from [ExAC](http://exac.broadinstitute.org/) |
| 101 - ExAC_AF_Adj | Adjusted Global Allele Frequency from ExAC |
| 102 - ExAC_AF_AFR | African/African American Allele Frequency from ExAC |
| 103 - ExAC_AF_AMR | American Allele Frequency from ExAC |
| 104 - ExAC_AF_EAS | East Asian Allele Frequency from ExAC |
| 105 - ExAC_AF_FIN | Finnish Allele Frequency from ExAC |
| 106 - ExAC_AF_NFE | Non-Finnish European Allele Frequency from ExAC |
| 107 - ExAC_AF_OTH | Other Allele Frequency from ExAC |
| 108 - ExAC_AF_SAS | South Asian Allele Frequency from ExAC |
| 109 - GENE_PHENO | Indicates if gene that the variant maps to is associated with a phenotype, disease or trait (0, 1, or null) |
| 110 - FILTER | Copied from input VCF |
| 111 - CONTEXT | The reference allele per VCF specs, and its 5 flanking base pairs |
| 112 - src_vcf_id | GDC UUID for the input VCF file |
| 113 - tumor_bam_uuid | GDC UUID for the tumor bam file |
| 114 - normal_bam_uuid | GDC UUID for the normal bam file |
| 115 - case_id | GDC UUID for the case |
| 116 - GDC_FILTER | GDC filters applied universally across all MAFs |
| 117 - COSMIC | Overlapping COSMIC variants |
| 118 - MC3_Overlap| Indicates whether this region overlaps with an MC3 variant for the same sample pair |
| 119 - GDC_Validation_Status | GDC implementation of validation checks |
| 120 - GDC_Valid_Somatic | True or False |
| 121 - vcf_region | Colon separated string containing the CHROM, POS, ID, REF, and ALT columns from the VCF file (e.g., chrZ:20:rs1234:A:T)|
| 122 - vcf_info | INFO column from VCF (not in somatic MAF)|
| 123 - vcf_format | FORMAT column from VCF (not in somatic MAF)|
| 124 - vcf_tumor_gt | Tumor sample genotype column from VCF (not in somatic MAF)|
| 125 - vcf_normal_gt | Normal sample genotype column from VCF (not in somatic MAF)|  
