# Submission Best Practices

Because of the data types and relationships included in the GDC, data submission can become a complex procedure. The purpose of this section is to present guidelines that will aid in the incorporation and harmonization of submitters' data. Please contact the GDC Help Desk at __<support@nci-gdc.datacommons.io>__ if you have any questions or concerns regarding a submission project.

## Date Obfuscation

The GDC is committed to providing accurate and useful information as well as protecting the privacy of  patients if necessary. Following this, the GDC accepts time intervals that were transformed to remove information that could identify an individual but preserve clinically useful timelines. The GDC recommends following a series of HIPAA regulations regarding the reporting of age-related information, which can be downloaded [here](BCRInformatics-DateObfuscator-291116-1100-2.pdf) as a PDF.

### General Guidelines

Actual calendar dates are not reported in GDC clinical fields but the lengths of time between events are preserved. Time points are reported based on the number of days since the patient's initial diagnosis. Events that occurred after the initial diagnosis are reported as positive and events that occurred before are reported as negative. Dates are not automatically obfuscated by the GDC validation system and submitters are required to make these changes in their clinical data. This affects these fields: `days_to_birth`, `days_to_death`, `days_to_last_follow_up`, `days_to_last_known_disease_status`, `days_to_recurrence`, `days_to_treatment`

>__Note:__ The day-based fields take leap years into account.  

### Patients Older than 90 Years and Clinical Events

Because of the low population number within the demographic of patients over 90 years old, it becomes more likely that patients can potentially be identified by a combination of their advanced age and publicly available clinical data. Because of this, patients over 90 years old are reported as exactly 90 years or 32,872 days old.

Following this, clinical events that occur over 32,872 days are also capped at 32,872 days. When timelines are capped, the priority should be to shorten the post-diagnosis values to preserve the accuracy of the age of the patient (except for patients who were diagnosed at over 90 years old). Values such as `days_to_death` and `days_to_recurrence` should be compressed before `days_to_birth` is compressed.

### Examples Timelines

__Example 1:__ An 88 year old patient is diagnosed with cancer and dies 13 years later.  The `days_to_birth` value is less than 32,872 days, so it can be accurately reported.  However, between the initial diagnosis and death, the patient turned 90 years old. Since 32,872 is the maximum, `days_to_death` would be calculated as 32872 - 32142 = 730.

>__Dates__

>* _Date of Birth:_ 01-01-1900
>* _Date of Initial Diagnosis:_ 01-01-1988
>* _Date of Death:_ 01-01-2001

>__Actual-Values__

>* _days_to_birth:_ -32142
>* _days_to_death:_ 4748

>__Obfuscated-Values__

>* _days_to_birth:_ -32142
>* _days_to_death:_ 730


__Example 2:__ A 98 year old patient is diagnosed with cancer and dies three years later.  Because `days_to_X` values are counted from initial diagnosis, days will be at their maximum value of 32,872 upon initial diagnosis. This will compress the later dates and reduce `days_to_birth` to -32,872 and `days_to_death` to zero.  

>__Dates__

>* _Date of Birth:_ 01-01-1900
>* _Date of Initial Diagnosis:_ 01-01-1998
>* _Date of Death:_ 01-01-2001

>__Actual-Values__

>* _days_to_birth:_ -35794
>* _days_to_death:_ 1095

>__Obfuscated-Values__

>* _days_to_birth:_ -32872
>* _days_to_death:_ 0

## Submitting Complex Data Model Relationships

The GDC Data Model includes relationships in which more than one entity of one type can be associated with one entity of another type. For example, more than one `read_group` entity can be associated with a `submitted_aligned_reads` entity. JSON-formatted files, in which a list object can be used, are well-suited to represent this type of relationship. Tab-delimited (TSV) files require additional syntax to demonstrate these relationships. For example, associating a `submitted_aligned_reads` entity to three read groups would require three `read_groups.submitter_id` columns, each with the `#` symbol and a number appended to them. See the two files below:

```TSV
type    submitter_id    data_category   data_format data_type   experimental_strategy   file_name   file_size   md5sum  read_groups.submitter_id#1 read_groups.submitter_id#2  read_groups.submitter_id#3
submitted_aligned_reads Alignment.bam  Raw Sequencing Data BAM Aligned Reads   WGS test_alignment.bam    123456789  aa6e82d11ccd8452f813a15a6d84faf1    READ_GROUP_1  READ_GROUP_2  READ_GROUP_3  

```
```JSON
{
    "type": "submitted_aligned_reads",
    "submitter_id": "Alignment.bam",
    "data_category": "Raw Sequencing Data",
    "data_format": "BAM",
    "data_type": "Aligned Reads",
    "experimental_strategy": "WGS",
    "file_name": "test_alignment.bam",
    "file_size": 123456789,
    "md5sum": "aa6e82d11ccd8452f813a15a6d84faf1",
    "read_groups": [
        {"submitter_id": "READ_GROUP_1"},
        {"submitter_id": "READ_GROUP_2"},
        {"submitter_id": "READ_GROUP_3"}
    ]
}
```

### Read groups

#### Submitting Read Group Names

The `read_group` entity requires a `read_group_name` field for submission.  If the `read_group` entity is associated with a BAM file, the submitter should use the `@RG` ID present in the BAM header as the `read_group_name`. This is important for the harmonization process and will reduce the possibility of errors.  

#### Multiple FASTQs from One Read Group

To align reads according to their direction and pair, the GDC requires that unaligned forward and reverse reads are submitted as "submitted_unaligned_reads."  When more than one FASTQ exists for a read group direction, the GDC requires that the FASTQ files are concatenated for each direction. In other words, each paired-end read group should be associated with exactly two FASTQ files (submitted_unaligned_reads).  

#### Minimal and Recommended Read Group Information

In addition to the required properties on `read_group` we also recommend submitting `flow_cell_barcode`, `lane_number` and `multiplex_barcode`.  This information can be used by our bioinformatics team and data downloaders to construct a `Platform Unit` (`PU`), which is a universally unique identifier that can be used to model various sequencing technical artifacts.  More information can be found in the [SAM specification PDF](https://github.com/samtools/hts-specs/blob/master/SAMv1.pdf).

For projects with library strategies of targeted sequencing or WXS we also require information on the target capture protocol included on `target_capture_kit`.

If this information is not provided it may cause a delay in the processing of submitted data.

Additional read group information will benefit data users.  Such information can be used by bioinformatics pipelines and will aid understanding and mitigation of batch effects.  If available, you should also provide as many of the remaining read group properties as possible.

## Submission File Quality Control

The GDC harmonization pipelines include multiple quality control steps to ensure the integrity of harmonized files and the efficient use of computational resources. For fastest turnaround of data processing we recommend that submitters perform basic QC of their data files prior to upload to identify any underlying data quality issues. This may include such tests as verifying expected genome coverage levels and sequencing quality.

## Target Capture Kit Q and A

1. What is a Target Capture Kit?  
Target capture kits contain reagents designed to enrich for and thus increase sequencing depth in certain genomic regions before library preparation. Two of the major methods to enrich genomic regions are by Hybridization and by PCR (Amplicon).

2. Why do we need Target Capture Kit information?  
Target region information is important for DNA-Seq variant calling and filtering, and essential for Copy Number Alternation and other analyses. This information is only needed for the Experimental Strategies of WXS or Targeted Sequencing.

3. How do submitters provide this information?  
There are 3 steps  
    * Step 1. The submitter should contact GDC User Service about any new Target Capture Kits that do not already exist in the GDC Dictionary. The GDC Bioinformatics and User Services teams will work together with the submitter to create a meaningful name for the kit and import this name and Target Region Bed File into the GDC.
    * Step 2. The submitter can then select one and only one GDC Target Capture Kit for each read group during molecular data submission.
    * Step 3. The submitter should also select the appropriate `library_strategy` and `library_selection` on the read_group entity.

4. What is a Target Region Bed File?  
A Target Region Bed File is tab-delimited file describing the kit target region in [bed format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1). The first 3 columns of such files are chrom, chromStart, and chromEnd.
Note that by definition, bed files are 0-based or "left-open, right-closed", which means bed interval "chr1    10    20" only contains 10 bases on chr1, from the 11th to the 20th.
In addition, submitters should also let GDC know the genome build (hg18, hg19 or GRCh38) of their bed files.

5. Is a Target Capture Kit uniquely defined by its Target Region Bed File?  
Not necessarily. Sometimes, users or manufactures may want to augment an existing kit with additional probes, in order to capture more regions or simply improve the quality of existing regions. In the latter case, the bed file stays the same, but it is now a different Target Capture Kit and should be registered separately as described in Step 3 above.

## Specifying Tumor Normal Pairs for Analysis

It is critical for many cancer bioinformatics pipelines to specify which normal sample to use to factor out germline variation.  In particular, this is a necessary specification for all tumor normal paired variant calling pipelines.  The following details describe how the GDC determines which normal sample to use for variant calling.

*  Every tumor aliquot will be used for variant calling.  For example, if 10 WXS tumor aliquots are submitted, the GDC will produce 10 alignments and 10 VCFs for each variant calling pipeline.
*  If there is only one normal we will use that normal for variant calling
*  If there are multiple normals of the same experimental_strategy for a case:
    *  Users can specify which normal to use by specifying on the aliquot.  To do so one of the following should be set to `TRUE` for the specified experimental strategy: `selected_normal_low_pass_wgs`, `selected_normal_targeted_sequencing`, `selected_normal_wgs`, or `selected_normal_wxs`.
    *  Or if no normal is specified the GDC will select the best normal for that patient based on the following criteria.  This same logic will also be used if multiple normal are selected.
        * If a case has blood cancer we will use sample type in the following priority order:

            Blood Derived Normal > Bone Marrow Normal > Mononuclear Cells from Bone Marrow Normal > Fibroblasts from Bone Marrow Normal > Lymphoid Normal > Buccal Cell Normal > Solid Tissue Normal > EBV Immortalized Normal

        * If a case does not have blood cancer we will use sample type in the following priority order:

            Solid Tissue Normal > Buccal Cell Normal > Lymphoid Normal > Fibroblasts from Bone Marrow Normal > Mononuclear Cells from Bone Marrow Normal > Bone Marrow Normal > Blood Derived Normal > EBV Immortalized Normal

        * If there are still ties, we will choose the aliquot submitted first.
* If there are no normals.
    * The GDC will not run tumor only variant calling pipeline by default.  The submitter must specify one of the following properties as TRUE: `no_matched_normal_low_pass_wgs`, `no_matched_normal_targeted_sequencing`, `no_matched_normal_wgs`, `no_matched_normal_wxs`.

Note that we will only run variant calling for a particular tumor aliquot per experimental strategy once. You must make sure that the appropriate normal control is uploaded to the GDC when Requesting Submission. Uploading a different normal sample later will not result in reanalysis by the GDC.

## Clinical Data Requirements

For the GDC to release a project there is a minimum number of clinical properties that are required.  Minimal cross-project GDC requirements include age, gender, and diagnosis information.  Other requirements may be added when the submitter is approved for submission to the GDC.

## miRNA Submission

The GDC requires that miRNA reads be trimmed before being uploaded to the GDC because miRNA datasets can have different trimming schemas. Uploading untrimmed miRNA reads can delay harmonization until the problem is resolved.  
