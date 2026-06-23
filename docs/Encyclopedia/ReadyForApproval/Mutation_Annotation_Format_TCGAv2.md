Mutation Annotation Format (MAF) - Legacy TCGA Specification
==============================================


*This definition was taken from the previously public wiki hosted by TCGA and reflects the MAF format
that was available during the active period of the TCGA project.*  




**Document Information**

The spec has been reverted to the June 26th version (version 20). Additional
changes are the removal of the "under construction" banner, changing all text to
black, and fixing a typo in the link to the MAF 2.2 specification.

**Specification for Mutation Annotation Format**  
Version 2.4.1  
June 20, 2014

**Contents**

-   1 Current version changes

-   2 About MAF specifications

    -   2.1 Definition of open access MAF
        data

    -   2.2 Somatic MAF vs. Protected
        MAF

-   3 MAF file fields

    -   3.1 Table 1 - File column
        headers

-   4 MAF file checks

-   5 MAF naming convention

-   6 Previous specification
    versions

Current version changes
=======================

This current revision is **version 2.4.1** of the Mutation Annotation Format
(MAF) specification.

The following items in the specification were added or modified in version 2.4.1
from version 2.4:

-   Header for MAF file is "\#version 2.4.1"

-   "Somatic" and "None" are the only acceptable values for "Mutation_Status"
    for a somatic.MAF (named .somatic.maf). When Mutation_Status is None,
    Validation_Status must be Invalid.

-   Centers need to make sure that Mutations_Status "None" doesn't include
    germline mutation.

-   For a somatic MAF, following rules should be satisfied:  
    SOMATIC = (A AND (B OR C OR D)) OR (E AND F)  
    A: *Mutation_Status* == "Somatic"  
    B: *Validation_Status* == "Valid"  
    C. *Verification_Status* == "Verified"  
    D. *Variant_Classification* is not {Intron, 5'UTR, 3'UTR, 5'Flank, 3'Flank,
    IGR}, which implies that *Variant_Classification* can only be
    \\{Frame_Shift_Del, Frame_Shift_Ins, In_Frame_Del, In_Frame_Ins,
    Missense_Mutation, Nonsense_Mutation, Silent, Splice_Site,
    Translation_Start_Site, Nonstop_Mutation, RNA, Targeted_Region}.  
    E: *Mutations_status == "None"*  
    F: *Validation_status == "Invalid"*

-   Extra validation rules: If Validation_Status == Valid or Invalid, then
    Validation_Method != none (case insensitive).

About MAF specifications
========================

Mutation annotation files should be transferred to the DCC. Those files should
be formatted using the mutation annotation format (MAF) that is described below.
File naming convention is also
[below](#MutationAnnotationFormat(MAF)Specificat).

Following categories of somatic mutations are reported in MAF files:

-   Missense and nonsense

-   Splice site, defined as SNP within 2 bp of the splice junction

-   Silent mutations

-   Indels that overlap the coding region or splice site of a gene or the
    targeted region of a genetic element of interest.

-   Frameshift mutations

-   Mutations in regulatory regions

### Definition of open access MAF data

A large proportion of MAFs are submitted as discovery data and sites labeled as
somatic in these files overlap with known germline variants. In order to
minimize germline contamination in putative (unvalidated) somatic calls, certain
filtering criteria have been imposed. Based on current policy, open access MAF
data should:

-   **include** all validated somatic mutation calls

-   **include** all unvalidated somatic mutation calls that overlap with a
    coding region or splice site

-   **exclude** all other types of mutation calls (i.e., non-somatic calls
    (validated or not), unvalidated somatic calls that are not in coding region
    or splice sites, and dbSNP sites that are not annotated as somatic in dbSNP,
    COSMIC or OMIM)



### Somatic MAF vs. Protected MAF

Centers will submit to the DCC MAF archives that contain Somatic MAF
(named**.somatic.maf**) for open access data and an all-inclusive Protected MAF
(named**.protected.maf**) that does not filter any data out and represents the
original super-set of mutation calls. The files will be formatted using the
Mutation Annotation Format (MAF).

The following table lists some of the critical attributes of somatic and
protected MAF files and provides a comparison.

| Attribute       | Somatic MAF | Protected MAF |
| -----------     | ----------- | ------------- |
| **File naming** | Somatic MAFs should be named as**\*.somatic.maf**and cannot contain 'germ' or 'protected' in file name. | Protected MAFs should be named as**\*.protected.maf**and should not contain 'somatic' in the file name. |
| **Mutation category**  | Somatic MAFs can only contain entries where*Mutation_Status*is "Somatic". If any other value is assigned to the field, the archive will fail. Experimentally validated or unvalidated (see next row) somatic mutations can be included in the file. | There is no such restriction for protected MAF. The file should contain all mutation calls including those from which .somatic.maf is derived. |
| **Filtering criteria** | In order to minimize germline contamination, somatic MAFs can contain unvalidated somatic mutations only from coding regions and splice sites, which implies: | There are no such constraints for mutations in protected MAF. |
|                        | If *Validation_Status* **is**"Unknown",*V a riant_Classification* **cannot** be 3'UTR, 3'Flank, 5'UTR, 5'Flank, IGR, or Intron.*Variant_Classification*can only be \\{Frame_Shift_Del, Frame_Shift_Ins, In_Frame_Del, In_Frame_Ins, Missense_Mutation, Nonsense_Mutation, Silent, Splice_Site, Translation_Start_Site, Nonstop_Mutation, RNA, Targeted_Region, De_novo_Start_InFrame, De_novo_Start_OutOfFrame\\}.                                                  |                                                                                                                                                |
|                        | There is no such constraint for experimentally validated (*Validation_Status*is "Valid") somatic mutations.                                                                  |                                                                                                                                                |
|                        |                                                                                |                                                                                                                                                |
|   | dbSNP sites that are not annotated as somatic in dbSNP, COSMIC or OMIM must be removed from somatic MAFs.                                                                         |                                                                                                                                                |
| **Access level**       | These files are deployed as open access data.    | These files are deployed as protected data.                                                                                                    |

MAF file fields
===============

The format of a MAF file is tab-delimited columns. Those columns are described
in Table 1 and are required in every MAF file. The order of the columns will be
validated by the DCC. Column headers and values **are** case sensitive where
specified. Columns may allow null values (i.e.\_ blank cells) and/or have
enumerated values. **The validator looks for a header stating the version of the
specification to validate against (e.g. \#version 2.4). If not, validation
fails.** Any columns that come after the columns described in Table 1 are
optional. Optional columns are not validated by the DCC and can be in any order.



Table 1 - File column headers
-----------------------------



| **Index** | **MAF Column Header** | **Description of Values** | **Example** | **Case Sensitive** | **Null** | **Enumerated** |                                    
| --------- | --------------------- | ------------------------- | ----------- | ------------------ | -------- | -------------- |
| 1         | Hugo_Symbol           | HUGO symbol for the gene (HUGO symbols are *always* in all caps). If no gene exists within 3kb enter "Unknown".  |EGFR                                           | Yes                   | No                                                                                                                                  | Set or Unknown                                                                                                                                                                                                                          |                                      |   |   |   |   |   |   |   |   |
|              |                               | Source: <http://genenames.org>                                                                                                                                                                                                                                                                                         |                                                |                       |                                                                                                                                     |                                                                                                                                                                                                                                         |                                      |   |   |   |   |   |   |   |   |
| 2            | Entrez_Gene_Id                | Entrez gene ID (an integer). If no gene exists within 3kb enter "0".                                                                                                                                                                                                                                                   | 1956                                           | No                    | No                                                                                                                                  | Set                                                                                                                                                                                                                                     |                                      |   |   |   |   |   |   |   |   |
|              |                               | Source: <http://ncbi.nlm.nih.gov/sites/entrez?db=gene>                                                                                                                                                                                                                                                                 |                                                |                       |                                                                                                                                     |                                                                                                                                                                                                                                         |                                      |   |   |   |   |   |   |   |   |
| 3            | Center                        | Genome sequencing center reporting the variant. If multiple institutions report the same mutation separate list using semicolons. Non-GSC centers will be also supported if center name is an accepted center name.                                                                                                    | hgsc.bcm.edu;genome.wustl.edu                  | Yes                   | No                                                                                                                                  | Set                                                                                                                                                                                                                                     |                                      |   |   |   |   |   |   |   |   |
| 4            | NCBI_Build                    | Any TGCA accepted genome identifier. Can be string, integer or a float.                                                                                                                                                                                                                                               | hg18, hg19, GRCh37, GRCh37-lite, 36, 36.1, 37, | No                    | No                                                                                                                                  | Set and Enumerated.                                                                                                                                                                                                                     |                                      |   |   |   |   |   |   |   |   |
| 5            | Chromosome                    | Chromosome number without "chr" prefix that contains the gene.                                                                                                                                                                                                                                                         | X, Y, M, 1, 2, etc.                            | Yes                   | No                                                                                                                                  | Set                                                                                                                                                                                                                                     |                                      |   |   |   |   |   |   |   |   |
| 6            | Start_Position                | Lowest numeric position of the reported variant on the genomic reference sequence. Mutation start coordinate (1-based coordinate system).                                                                                                                                                                              | 999                                            | No                    | No                                                                                                                                  | Set                                                                                                                                                                                                                                     |                                      |   |   |   |   |   |   |   |   |
| 7            | End_Position                  | Highest numeric genomic position of the reported variant on the genomic reference sequence. Mutation end coordinate (inclusive, 1-based coordinate system).                                                                                                                                                            | 1000                                           | No                    | No                                                                                                                                  | Set                                                                                                                                                                                                                                     |                                      |   |   |   |   |   |   |   |   |
| 8            | Strand                        | Genomic strand of the reported allele. Variants should always be reported on the positive genomic strand. (Currently, only the positive strand is an accepted value).                                                                                                                                                  | \+                                             | No                    | No                                                                                                                                  | \+                                                                                                                                                                                                                                      |                                      |   |   |   |   |   |   |   |   |
| 9            | Variant_Classification        | Translational effect of variant allele.                                                                                                                                                                                                                                                                                | Missense_Mutation                              | Yes                   | No                                                                                                                                  | Frame_Shift_Del, Frame_Shift_Ins, In_Frame_Del, In_Frame_Ins, Missense_Mutation, Nonsense_Mutation, Silent, Splice_Site, Translation_Start_Site, Nonstop_Mutation, 3'UTR, 3'Flank, 5'UTR, 5'Flank, IGR *(See Notes Section #1)* , Intron, RNA, Targeted_Region |                                      |   |   |   |   |   |   |   |   |
| 10           | Variant_Type                  | Type of mutation. TNP (tri-nucleotide polymorphism) is analogous to DNP but for 3 consecutive nucleotides. ONP (oligo-nucleotide polymorphism) is analogous to TNP but for consecutive runs of 4 or more.                                                                                                              | INS                                            | Yes                   | No                                                                                                                                  | SNP, DNP, TNP, ONP, INS, DEL, or Consolidated *(See Notes Section #2)* )                                                                                                                                                                                         |                                      |   |   |   |   |   |   |   |   |
| 11           | Reference_Allele              | The plus strand reference allele at this position. Include the sequence deleted for a deletion, or "-" for an insertion.                                                                                                                                                                                               | A                                              | Yes                   | No                                                                                                                                  | A,C,G,T and/or -                                                                                                                                                                                                                        |                                      |   |   |   |   |   |   |   |   |
| 12           | Tumor_Seq_Allele1             | Primary data genotype. Tumor sequencing (discovery) allele 1. " -" for a deletion represent a variant. "-" for an insertion represents wild-type allele. Novel inserted sequence for insertion should not include flanking reference bases.                                                                            | C                                              | Yes                   | No                                                                                                                                  | A,C,G,T and/or -                                                                                                                                                                                                                        |                                      |   |   |   |   |   |   |   |   |
| 13           | Tumor_Seq_Allele2             | Primary data genotype. Tumor sequencing (discovery) allele 2. " -" for a deletion represents a variant. "-" for an insertion represents wild-type allele. Novel inserted sequence for insertion should not include flanking reference bases.                                                                           | G                                              | Yes                   | No                                                                                                                                  | A,C,G,T and/or -                                                                                                                                                                                                                        |                                      |   |   |   |   |   |   |   |   |
| 14           | dbSNP_RS                      | Latest dbSNP rs ID (dbSNP_ID) or "novel" if there is no dbSNP record. source: <http://ncbi.nlm.nih.gov/projects/SNP/>                                                                                                                                                                                                  | rs12345                                        | Yes                   | Yes                                                                                                                                 | Set or "novel"                                                                                                                                                                                                                          |                                      |   |   |   |   |   |   |   |   |
| 15           | dbSNP_Val_Status              | dbSNP validation status. Semicolon- separated list of validation statuses.                                                                                                                                                                                                                                             | by2Hit2Allele;byCluster                        | No                    | Yes                                                                                                                                 | by1000genomes;by2Hit2Allele; byCluster; byFrequency; byHapMap; byOtherPop; bySubmitter; alternate_allele *(See Notes Section #3)* **Note that "none" will no longer be an acceptable value.**                                                                  |                                      |   |   |   |   |   |   |   |   |
| 16           | Tumor_Sample_Barcode          | BCR aliquot barcode for the tumor sample including the two additional fields indicating plate and well position. i.e. TCGA-SiteID-PatientID-SampleID-PortionID-PlateID-CenterID. The full TCGA Aliquot ID.                                                                                                             | TCGA-02-0021-01A-01D-0002-04                   | Yes                   | No                                                                                                                                  | Set                                                                                                                                                                                                                                     |                                      |   |   |   |   |   |   |   |   |
| 17           | Matched_Norm_Sample_Barcode   | BCR aliquot barcode for the matched normal sample including the two additional fields indicating plate and well position. i.e. TCGA-SiteID-PatientID-SampleID-PortionID-PlateID-CenterID. The full TCGA Aliquot ID; e.g. TCGA-02-0021-10A-01D-0002-04 (compare portion ID '10A' normal sample, to '01A' tumor sample). | TCGA-02-0021-10A-01D-0002-04                   | Yes                   | No                                                                                                                                  | Set                                                                                                                                                                                                                                     |                                      |   |   |   |   |   |   |   |   |
| 18           | Match_Norm_Seq_Allele1        | Primary data. Matched normal sequencing allele 1. "-" for deletions; novel inserted sequence for INS not including flanking reference bases.                                                                                                                                                                           | T                                              | Yes                   | Yes                                                                                                                                 | A,C,G,T and/or -                                                                                                                                                                                                                        |                                      |   |   |   |   |   |   |   |   |
| 19           | Match_Norm_Seq_Allele2        | Primary data. Matched normal sequencing allele 2. "-" for deletions; novel inserted sequence for INS not including flanking reference bases.                                                                                                                                                                           | ACGT                                           | Yes                   | Yes                                                                                                                                 | A,C,G,T and/or -                                                                                                                                                                                                                        |                                      |   |   |   |   |   |   |   |   |
| 20           | Tumor_Validation_Allele1      | Secondary data from orthogonal technology. Tumor genotyping (validation) for allele 1. "-" for deletions; novel inserted sequence for INS not including flanking reference bases.                                                                                                                                      | \-                                             | Yes                   | Yes                                                                                                                                 | A,C,G,T and/or -                                                                                                                                                                                                                        |                                      |   |   |   |   |   |   |   |   |
| 21           | Tumor_Validation_Allele2      | Secondary data from orthogonal technology. Tumor genotyping (validation) for allele 2. "-" for deletions; novel inserted sequence for INS not including flanking reference bases.                                                                                                                                      | A                                              | Yes                   | Yes                                                                                                                                 | A,C,G,T and/or -                                                                                                                                                                                                                        |                                      |   |   |   |   |   |   |   |   |
| 22           | Match_Norm_Validation_Allele1 | Secondary data from orthogonal technology. Matched normal genotyping (validation) for allele 1. "-" for deletions; novel inserted sequence for INS not including flanking reference bases.                                                                                                                             | C                                              | Yes                   | Yes                                                                                                                                 | A,C,G,T and/or -                                                                                                                                                                                                                        |                                      |   |   |   |   |   |   |   |   |
| 23           | Match_Norm_Validation_Allele2 | Secondary data from orthogonal technology. Matched normal genotyping (validation) for allele 2. "-" for deletions; novel inserted sequence for INS not including flanking reference bases.                                                                                                                             | G                                              | Yes                   | Yes                                                                                                                                 | A,C,G,T and/or -                                                                                                                                                                                                                        |                                      |   |   |   |   |   |   |   |   |
| 24           | Verification_Status *(See Notes Section #4)*         | Second pass results from independent attempt using same methods as primary data source. Generally reserved for 3730 Sanger Sequencing.                                                                                                                                                                                 | Verified                                       | Yes                   | Yes                                                                                                                                 | Verified, Unknown                                                                                                                                                                                                                       |                                      |   |   |   |   |   |   |   |   |
| 25           | Validation_Status *(See Notes Section #5)*           | Second pass results from orthogonal technology.                                                                                                                                                                                                                                                                        | Valid                                          | Yes                   | No                                                                                                                                  | Untested, Inconclusive, Valid, Invaild                                                                                                                                                                                                                                |                                      |   |   |   |   |   |   |   |   |
| 26           | Mutation_Status               | Updated to reflect validation or verification status and to be in agreement with the [VCF VLS](https://wiki.nci.nih.gov/x/2gcYAw) field. The values allowed in this field are constrained by the value in the Validation_Status field.                                                                                 | Somatic                                        | Yes                   | No                                                                                                                                  | **Validation_Status values:** Untested, Inconslusive, Valid, Invalid - **Allowed Mutations_Status Values for Untested and Inconclusive:** *(See Notes Seciton #6)* None, Germline, Somatic, LOH, Post-transcriptional modification **Unknown Allowed Mutation_status Values for Valid:** *(See Notes Seciton #6)* Germline, Somatic, LOH, Post-transcriptional modification, Unknown - **Allowed Mutations_Status Values for Invalid:** *(See Notes Seciton #6)* none                                                                                                                                                                                                             |  |   |   |   |   |   |   |   |   |
 |                                                                                                                                                                                                                                                                                                                        |                                                |                       |                                                                                                                                     |                                                                                                                                                                                                                                         |                                      |   |   |   |   |   |   |   |   |
| 27           | Sequencing_Phase              | TCGA sequencing phase. Phase should change under any circumstance that the targets under consideration change.                                                                                                                                                                                                         | Phase_I                                        | No                    | Yes                                                                                                                                 | No                                                                                                                                                                                                                                      |                                      |   |   |   |   |   |   |   |   |
| 28           | Sequence_Source               | Molecular assay type used to produce the analytes used for sequencing. Allowed values are a subset of the [SRA 1.5](http://www.ncbi.nlm.nih.gov/viewvc/v1/trunk/sra/doc/SRA_1-5/) library_strategy field values. This subset matches those used at CGHub.                                                              | WGS;WXS                                        | Yes                   | No                                                                                                                                  | **Common TCGA values:** WGS, WGA, WXS, RNA-Seq, miRNA-Seq, Bisulfite-Seq, VALIDATION, Other **Other allowed values (per SRA 1.5)**  ncRNA-Seq, WCS, CLONE, POOLCLONE, AMPLICON, CLONEEND, FINISHING, ChIP-Seq, MNase-Seq, DNase-Hypersensitivity, EST, FL-cDNA, CTS, MRE-Seq, MeDIP-Seq, MBD-Seq, Tn-Seq, FAIRE-seq, SELEX, RIP-Seq, ChIA-PET
                                                                                                                                                                                                                       |                                      |   |   |   |   |   |   |   |   |
| 29           | Validation_Method             | The assay platforms used for the validation call. Examples: Sanger_PCR_WGA, Sanger_PCR_gDNA, 454_PCR_WGA, 454_PCR_gDNA; separate multiple entries using semicolons.                                                                                                                                                    | Sanger_PCR_WGA;Sanger_PCR_gDNA                 | No                    | **NO**. I**f Validation_Status = Untested then "none"** If Validation_Status = Valid or Invalid, then not "none" (case insensitive) | No                                                                                                                                                                                                                                      |                                      |   |   |   |   |   |   |   |   |
| 30           | Score                         | Not in use.                                                                                                                                                                                                                                                                                                            | NA                                             | No                    | Yes                                                                                                                                 | No                                                                                                                                                                                                                                      |                                      |   |   |   |   |   |   |   |   |
| 31           | BAM_File                      | Not in use.                                                                                                                                                                                                                                                                                                            | NA                                             | No                    | Yes                                                                                                                                 | No                                                                                                                                                                                                                                      |                                      |   |   |   |   |   |   |   |   |
| 32           | Sequencer                     | Instrument used to produce primary data. Separate multiple entries using semicolons.                                                                                                                                                                                                                                   | Illumina GAIIx;SOLID                           | Yes                   | No                                                                                                                                  | Illumina GAIIx, Illumina HiSeq, SOLID, 454, ABI 3730xl, Ion Torrent PGM, Ion Torrent Proton, PacBio RS, Illumina MiSeq, Illumina HiSeq 2500, 454 GS FLX Titanium, AB SOLiD 4 System                                                                                                                                                                                                                         |                                      |   |   |   |   |   |   |   |   |
| 33           | Tumor_Sample_UUID             | BCR aliquot UUID for tumor sample                                                                                                                                                                                                                                                                                      | 550e8400-e29b-41d4-a716-446655440000           | Yes                   | No                                                                                                                                  |                                                                                                                                                                                                                                         |                                      |   |   |   |   |   |   |   |   |
| 34           | Matched_Norm_Sample_UUID      | BCR aliquot UUID for matched normal                                                                                                                                                                                                                                                                                    | 567e8487-e29b-32d4-a716-446655443246           | Yes                   | No                                                       |                                                                                                                                                  

**Notes**<br/>
*1 Intergenic Region.* <br/>
*2 Consolidationd is used to indicate a site that was initially reported as a variant but subsequently removed from further analysis because it was consolidated into a new variant. For example, a SNP variant incorporated into a TNP variant.* <br/>
*3 Used when the discovered varieant differs from that of dbSNP.*<br/>
*4 These MAF headers describe the technology that was used to confirm a mutation, whether the same technology ("verification") or a different technology ("validation") is used to prove that a variant is germline or a somatic mutation.*<br/>
*5 These MAF headers describe the technology that was used toconfirm a mutation, whether the same technology (verification) or a different technology (validation) is used to prove that a variant is germline or a somatic mutation.*<br/>
*6 Explanation of some Validation Status-Mutation Status combinations.*<br/>

| Validation Status | Mutation Status | Explanation |
| ------------------ | --------------- | ----------- |
| Valid              | Unknown         | a valid variant with unknown somatic status due to lack of data from matched normal tissue. |
| Invalid            | None            | validation attempted, tumor and normal are homozygous reference (formerly described as Wildtype) |
| Inconclusive       | Unknown         | validation failed, neither the genotype nor its somatic status is certain due to lack of data from matched normal tissue |
| Inconclusive       | None            | validation failed, tumor genotype appears to be homozygous reference |

    Important Criteria

    **Index column indicates the order in which the columns are expected**. **All
    headers are case sensitive.** The Case Sensitive column specifies which values
    are case sensitive. The Null column indicates which MAF columns are allowed to
    have null values. The Enumerated column indicates which MAF columns have
    specified values: an Enumerated value of "No" indicates that there are no
    specified values for that column; other values indicate the specific values
    listed allowed; a value of "Set" indicates that the MAF column values come from
    a specified set of known values (*e.g.*HUGO gene symbols).


MAF file checks
===============

The DCC Archive Validator checks the integrity of a MAF file. Validation will
fail if any of the below are not true for a MAF file:

1.  Column header text (including case) and order must match specification
    (Table 1) exactly

2.  Values under column headers listed in the specification (Table 1) as not
    null must have values

3.  Values that are specified in Table 1 as Case Sensitive must be.

4.  If column headers are listed in the specification as having *enumerated*
    values (*i.e.* a "Yes" in the "Enumerated" column), then the values under
    those column must come from the enumerated values listed under "Enumerated".

5.  If column headers are listed in the specification as having *set* values
    (*i.e.* a "Set" in the "Enumerated" column), then the values under those
    column must come from the enumerated values of that domain (*e.g.* HUGO gene
    symbols).

6.  All Allele-based columns must contain- (deletion), or a string composed of
    the following capitalized letters: A, T, G, C.

7.  IfValidation_Status== "Untested"
    thenTumor_Validation_Allele1,Tumor_Validation_Allele2,Match_Norm_Validation_Allele1,Match_Norm_Validation_Allele2can
    be null (depending onValidation_Status).

    1.  IfValidation_Status== "Inconclusive"
        thenTumor_Validation_Allele1,Tumor_Validation_Allele2,Match_Norm_Validation_Allele1,Match_Norm_Validation_Allele2can
        be null (depending onValidation_Status)**.**

8.  If Validation_Status == Valid, then Validated_Tumor_Allele1 and
    Validated_Tumor_Allele2must be populated (one of A, C, G, T, and -)

    1.  If Validation_Status == "Valid" then Tumor_Validation_Allele1,
        Tumor_Validation_Allele2, Match_Norm_Validation_Allele1,
        Match_Norm_Validation_Allele2 cannot be null

    2.  IfValidation_Status== "Invalid"
        thenTumor_Validation_Allele1,Tumor_Validation_Allele2,Match_Norm_Validation_Allele1,Match_Norm_Validation_Allele2cannot
        be null AND Tumor_Validation_Allelle1 ==
        Match_Norm_Validation_Allele1AND Tumor_Validation_Allelle2 ==
        Match_Norm_Validation_Allele2 (Added as a replacement for 8a as a
        result of breakdown)

9.  Check allele values against Mutation_Status:  
    Check allele values against Validation_status:

    1.  If Mutation_Status == "Germline" and Validation_Status == "Valid", then
        Tumor_Validation_Allele1 == Match_Norm_Validation_Allele1 and
        Tumor_Validation_Allele2 == Match_Norm_Validation_Allele2.

    2.  If Mutation_Status == "Somatic" and Validation_Status == "Valid", then
        Match_Norm_Validation_Allele1 == Match_Norm_Validation_Allele2 ==
        Reference_Allele and (Tumor_Validation_Allele1 or
        Tumor_Validation_Allele2) != Reference_Allele

    3.  If Mutation_Status == "LOH" and Validation_Status=="Valid", then
        Tumor_Validation_Allele1 == Tumor_Validation_Allele2 and
        Match_Norm_Validation_Allele1 != Match_Norm_Validation_Allele2 and
        Tumor_Validation_Allele1 == (Match_Norm_Validation_Allele1 or
        Match_Norm_Validation_Allele2).

10. Check that Start_position \<= End_position

11. Check for the Start_position and End_position against Variant_Type:

    1.  If Variant_Type == "INS", then (End_position - Start_position + 1 ==
        length (Reference_Allele) or End_position - Start_position == 1) and
        length(Reference_Allele) \<= length(Tumor_Seq_Allele1 and
        Tumor_Seq_Allele2)

    2.  If Variant_Type == "DEL", then End_position - Start_position + 1 ==
        length (Reference_Allele), then length(Reference_Allele) \>=
        length(Tumor_Seq_Allele1 and Tumor_Seq_Allele2)

    3.  If Variant_Type == "SNP", then length(Reference_Allele and
        Tumor_Seq_Allele1 and Tumor_Seq_Allele2) == 1 and (Reference_Allele and
        Tumor_Seq_Allele1 and Tumor_Seq_Allele2) != "-"

    4.  If Variant_Type == "DNP", then length(Reference_Allele and
        Tumor_Seq_Allele1 and Tumor_Seq_Allele2) == 2 and (Reference_Allele and
        Tumor_Seq_Allele1 and Tumor_Seq_Allele2) !contain "-"

    5.  If Variant_Type == "TNP", then length(Reference_Allele and
        Tumor_Seq_Allele1 and Tumor_Seq_Allele2) == 3 and (Reference_Allele and
        Tumor_Seq_Allele1 and Tumor_Seq_Allele2) !contain "-"

    6.  If Variant_Type == "ONP", then length(Reference_Allele) ==
        length(Tumor_Seq_Allele1) == length(Tumor_Seq_Allele2) \> 3 and
        (Reference_Allele and Tumor_Seq_Allele1 and Tumor_Seq_Allele2) !contain
        "-"

12. Validation for UUID-based files:

    1.  Column \#33 must be Tumor_Sample_UUID containing UUID of the BCR aliquot
        for tumor sample

    2.  Column \#34 must be Matched_Norm_Sample_UUID containing UUID of the BCR
        aliquot for matched normal sample

    3.  Metadata represented by Tumor_Sample_Barcode and
        Matched_Norm_Sample_Barcode should correspond to the UUIDs assigned to
        Tumor_Sample_UUID and Matched_Norm_Sample_UUID respectively

13. If Validation_Status == "Valid" or "Invalid", then Validation_Method !=
    "none" (case insensitive) .

MAF naming convention
=====================

In archives uploaded to the DCC, the MAF file name should relate to the
containing archive name in the following way:

If the archive has the name

    \<domain\>_\<disease_abbrev\>.\<platform\>.Level_2.\<serial_index\>.\<revision\>.0.tar.gz

then a somatic MAF file with the archive should be named according to

    \<domain\>_\<disease_abbrev\>.\<platform\>.Level_2.\<serial_index\>[.\<optional_tag\>].somatic.maf

and a protected MAF with the archive should be named according to

    \<domain\>_\<disease_abbrev\>.\<platform\>.Level_2.\<serial_index\>[.\<optional_tag\>].protected.maf

The \<optional_tag\> may consist of alphanumeric characters, dash, and
underscore; no spaces or periods; or it may be left out altogether. The purpose
of the optional tag is to impart some brief annotation.

*Example*

For the archive

    genome.wustl.edu_OV.IlluminaGA_DNASeq.Level_2.7.6.0.tar.gz

the following are examples of valid maf names

    genome.wustl.edu_OV.IlluminaGA_DNASeq.Level_2.7.somatic.maf
    genome.wustl.edu_OV.IlluminaGA_DNASeq.Level_2.7.protected.maf
