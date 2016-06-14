# GDC MAF Format

## Introduction

Mutation Annotation Format (MAF) is a tab-delimited text file that lists mutations.  The GDC produces two types of MAF files, Protected and Somatic (or Public) MAFs, for each variant calling pipeline in each GDC project. These MAFs are derived from the GDC annotated VCF files.

Annotated VCF files often have variants reported on multiple transcripts whereas the protected MAF (\*protected.maf) only reports the most critically affected one.  Somatic MAFs (\*somatic.maf) are further processed to remove low quality and potential germline variants. In addition, for tumor samples that contain variants from multiple combinations of tumor normal aliquot pairs, only one pair is selected in the Somatic MAF. Somatic MAFs are open access and can be freely distributed.

The GDC MAF file format follows the standard format of the <a href="https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification">TCGA Mutation Annotation Format</a>, but with additional columns.  


## Protected MAF File Structure

Listed below are the columns in a protected MAF. Most columns are explained in either <a href="https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification">TCGA MAF specification</a> or <a href="http://www.ensembl.org/info/docs/tools/vep/vep_formats.html">VEP documentation</a>.

### Important updates from TCGA v2.4

1. Column #4 “NCBI_Build” is GRCh38 by default
2. Column #60 VEP name "STRAND" is changed to "TRANSCRIPT_STRAND" to avoid confusion with Column#8 "Strand”
3. Column #15 “dbSNP_Val_Status" is not implemented
4. Column #26 “Mutation_Status”: multiallelic calls from MuSE are labeled in this column as  MuSEMulti
5. Column #32 “Sequencer” includes sequencers used. If different sequencers were used to generate normal and tumor data, normal sequencer is listed first.
6. Column #113-116 “vcf_info”, “vcf_format”, “vcf_tumor_gt” and “vcf_normal_gt” are the corresponding columns from VCF files. Including them facilitates parsing specific variant information.
7. Column #109 "FILTER" entries are separated by both commas and semi-colons.
8. Column #117 “GDC_Validation_Status”: GDC also collects TCGA validation sequences.  It compares these with variants derived from Next-Generation Sequencing data from the same sample and populates the comparison result in "GDC_Validation_Status".
    * "Valid", if Alternative allele(s) in the tumor validation sequence is(are) the same as GDC variant call
    * "Invalid", if none of the alternative allele(s) in the tumor validation sequence is the same as GDC variant call
    * "Inconclusive" if two alternative allele exists, and one matches while the other does not.
    * "Unknown" if no validation sequence exists
9. Column #118 GDC_Valid_Somatic is TRUE if GDC_Validation_Status is "Valid" and the variant is "Somatic" in validation calls.  It is FALSE if these criteria are not met.

### List of Columns
1. Hugo_Symbol
2. Entrez_Gene_Id
3. Center
4. NCBI_Build
5. Chromosome
6. Start_Position
7. End_Position
8. Strand
9. Variant_Classification
10. Variant_Type
11. Reference_Allele
12. Tumor_Seq_Allele1
13. Tumor_Seq_Allele2
14. dbSNP_RS
15. dbSNP_Val_Status
16. Tumor_Sample_Barcode
17. Matched_Norm_Sample_Barcode
18. Match_Norm_Seq_Allele1
19. Match_Norm_Seq_Allele2
20. Tumor_Validation_Allele1
21. Tumor_Validation_Allele2
22. Match_Norm_Validation_Allele1
23. Match_Norm_Validation_Allele2
24. Verification_Status
25. Validation_Status
26. Mutation_Status
27. Sequencing_Phase
28. Sequence_Source
29. Validation_Method
30. Score
31. BAM_File
32. Sequencer
33. Tumor_Sample_UUID
34. Matched_Norm_Sample_UUID
35. HGVSc
36. HGVSp
37. HGVSp_Short
38. Transcript_ID
39. Exon_Number
40. t_depth
41. t_ref_count
42. t_alt_count
43. n_depth
44. n_ref_count
45. n_alt_count
46. all_effects
47. Allele
48. Gene
49. Feature
50. Feature_type
51. Consequence
52. cDNA_position
53. CDS_position
54. Protein_position
55. Amino_acids
56. Codons
57. Existing_variation
58. ALLELE_NUM
59. DISTANCE
60. TRANSCRIPT_STRAND
61. SYMBOL
62. SYMBOL_SOURCE
63. HGNC_ID
64. BIOTYPE
65. CANONICAL
66. CCDS
67. ENSP
68. SWISSPROT
69. TREMBL
70. UNIPARC
71. RefSeq
72. SIFT
73. PolyPhen
74. EXON
75. INTRON
76. DOMAINS
77. GMAF
78. AFR_MAF
79. AMR_MAF
80. ASN_MAF
81. EAS_MAF
82. EUR_MAF
83. SAS_MAF
84. AA_MAF
85. EA_MAF
86. CLIN_SIG
87. SOMATIC
88. PUBMED
89. MOTIF_NAME
90. MOTIF_POS
91. HIGH_INF_POS
92. MOTIF_SCORE_CHANGE
93. IMPACT
94. PICK
95. VARIANT_CLASS
96. TSL
97. HGVS_OFFSET
98. PHENO
99. MINIMISED
100. ExAC_AF
101. ExAC_AF_AFR
102. ExAC_AF_AMR
103. ExAC_AF_EAS
104. ExAC_AF_FIN
105. ExAC_AF_NFE
106. ExAC_AF_OTH
107. ExAC_AF_SAS
108. GENE_PHENO
109. FILTER
110. src_vcf_id
111. tumor_bam_uuid
112. normal_bam_uuid
113. vcf_info
114. vcf_format
115. vcf_tumor_gt
116. vcf_normal_gt
117. GDC_Validation_Status
118. GDC_Valid_Somatic

## Somatic MAF file structure

This is very similar to the public MAF format, but with some important additional filters.  The process for modifying a Protected MAF into a Somatic MAF are the following:

1. Aliquot Selection: only one tumor-normal pair are selected for each tumor sample based on the plate number, sample_type, analyte_type and other features extracted from tumor TCGA aliquot barcode.
2. Low quality variants filtering: currently only FILTER="PASS" are kept.
Also removed: n_depth < 8
3. Germline Masking
    * Remove variants with Mutation_Status does not equal "Somatic"
    * Keep variants if any of the follow three criteria applies, but dbSNP sites that are not annotated as "Somatic" in VEP are also removed
        * GDC_Valid_Somatic is TRUE
        * VEP "Somatic" field is not empty (currently not implemented)
        * The variant affects transcript sequences.  Variant_Classification cannot be 3'UTR, 3'Flank, 5'UTR, 5'Flank, IGR, or Intron. Variant_Classification can only be Frame_Shift_Del, Frame_Shift_Ins, In_Frame_Del, In_Frame_Ins, Missense_Mutation, Nonsense_Mutation, Silent, Splice_Site, Translation_Start_Site, Nonstop_Mutation, RNA, Targeted_Region, De_novo_Start_InFrame, or De_novo_Start_OutOfFrame.
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
