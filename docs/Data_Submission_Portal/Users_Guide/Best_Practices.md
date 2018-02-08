# Submission Best Practices

Because of the data types and relationships included in the GDC, data submission can become a complex procedure. The purpose of this section is to present guidelines that will aid in the incorporation and harmonization of submitters' data. Please contact the GDC Help Desk at __support@nci-gdc.datacommons.io__ if you have any questions or concerns regarding a submission project.

## Date Obfuscation

The GDC is committed to providing accurate and useful information as well as protecting the privacy of  patients if necessary. Following this, the GDC accepts time intervals that were transformed to remove information that could identify an individual but preserve clinically useful timelines. The GDC recommends following a series of HIPAA regulations regarding the reporting of age-related information, which can be downloaded [here](BCRInformatics-DateObfuscator-291116-1100-2.pdf) as a PDF.

### General Guidelines

Actual calendar dates are not reported in GDC clinical fields but the lengths of time between events are preserved. Time points are reported based on the number of days since the patient's initial diagnosis. Events that occurred after the initial diagnosis are reported as positive and events that occurred before are reported as negative. Dates are not automatically obfuscated by the GDC validation system and submitters are required to make these changes in their clinical data.

| Affected Fields |
| --- |
| `days_to_birth` |
| `days_to_death` |
| `days_to_last_follow_up` |
| `days_to_last_known_disease_status` |
| `days_to_recurrence` |
| `days_to_treatment` |

### Patients Older than 90 Years

Because of the low population number within the demographic of patients over 90 years old, it becomes more likely that patients can potentially be identified by a combination of their advanced age and publicly available clinical data. Because of this, patients over 90 years old are reported as exactly 90 years or 32,872 days old.

__Note:__ The day-based fields take leap years into account.   

### Clinical Events After a Patient Turns 90 Years Old

Clinical events that occur over 32,872 days after an event also have the potential to reveal the age and identity of an individual over the age of 90. Following this, all timelines are capped at 32,872 days. When timelines are capped, the priority should be to shorten the post-diagnosis values to preserve the accuracy of the age of the patient (except for patients who were diagnosed at over 90 years old). Values such as `days_to_death` and `days_to_recurrence` should be compressed before `days_to_birth` is compressed.

### Examples Timelines

__Example 1:__ An 88 year old patient is diagnosed with cancer and dies 13 years later.  The `days_to_birth` value is less than 32,872 days, so it can be accurately reported.  However, between the initial diagnosis and death, the patient turned 90 years old. Since 32,872 is the maximum, `days_to_death` would be calculated as 32872 - 32142 = 730.

__Dates__

* _Date of Birth:_ 01-01-1900
* _Date of Initial Diagnosis:_ 01-01-1988
* _Date of Death:_ 01-01-2001

__Actual-Values__

* _days_to_birth:_ -32142
* _days_to_death:_ 4748

__Obfuscated-Values__

* _days_to_birth:_ -32142
* _days_to_death:_ 730


__Example 2:__ A 98 year old patient is diagnosed with cancer and dies three years later.  Because `days_to_X` values are counted from initial diagnosis, days will be at their maximum value of 32,872 upon initial diagnosis. This will compress the later dates and reduce `days_to_birth` to -32,872 and `days_to_death` to zero.  

__Dates__

* _Date of Birth:_ 01-01-1900
* _Date of Initial Diagnosis:_ 01-01-1998
* _Date of Death:_ 01-01-2001

__Actual-Values__

* _days_to_birth:_ -35794
* _days_to_death:_ 1095

__Obfuscated-Values__

* _days_to_birth:_ -32872
* _days_to_death:_ 0


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

## Read groups

### Submitting Read Group Names

The `read_group` entity requires a `read_group_name` field for submission.  If the `read_group` entity is associated with a BAM file, the submitter should use the `@RG` ID present in the BAM header as the `read_group_name`. This is important for the harmonization process and will reduce the possibility of errors.  

### Minimal Read Group Information

In addition to the required properties on `read_group` we also recommend submitting `flow_cell_barcode`, `lane_number` and `multiplex_barcode`.  This information can be used by our bioinformatics team and data downloaders to construct a `Platform Unit` (`PU`), which is a universally unique identifier that can be used to model various sequencing technical artifacts.  More information can be found in the SAM specification (https://github.com/samtools/hts-specs/blob/master/SAMv1.pdf).

For projects with library strategies of targeted sequencing or WXS we also require information on the target capture protocol or the following properties:

* `target_capture_kit_catalog_number`
* `target_capture_kit_name`
* `target_capture_kit_target_region`
* `target_capture_kit_vendor`
* `target_capture_kit_version`

If this information is not provided it may cause a delay in the processing of submitted data.

### Recommended Read Group Information

Additional read group information will benefit data users.  Such information can be used by bioinformatics pipelines and will aid understanding and mitigation of batch effects.  If available you should also provide as many of the remaining read group properties as possible.

## Submission File Quality Control

The GDC harmonization pipelines include multiple quality control steps to ensure the integrity of harmonized files and the efficient use of computational resources. For fastest turnaround of data processing we recommend that submitters perform basic QC of their data files prior to upload to identify any underlying data quality issues. This may include such tests as verifying expected genome coverage levels and sequencing quality.
