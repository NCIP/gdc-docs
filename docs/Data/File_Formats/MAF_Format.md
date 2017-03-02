# GDC MAF Format
- Version 1.0.0

## Introduction

Mutation Annotation Format (MAF) is a tab-delimited text file that describes mutations and are generated on a project-level.  The GDC produces MAF files at two permission levels: __protected__ and __somatic__ (or open-access). One MAF files is produced per variant calling pipeline per GDC project. MAFs are produced by aggregating the GDC annotated VCF files.

Annotated VCF files often have variants reported on multiple transcripts whereas the protected MAF file generated from the VCFs (\*protected.maf) only reports the most critically affected one. Somatic MAFs (\*somatic.maf) are further processed to remove low quality and potential germline variants. For tumor samples that contain variants from multiple combinations of tumor normal aliquot pairs, only one pair is selected in the Somatic MAF. Somatic MAFs are publicly available and can be freely distributed.

The GDC MAF file format is based on the [TCGA Mutation Annotation Format](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification) specifications, but some additional columns were included and some were removed.

## Protected MAF File Structure

Listed below are the columns in a protected MAF and their definitions.

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
