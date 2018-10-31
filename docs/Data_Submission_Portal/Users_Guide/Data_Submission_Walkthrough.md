# Data Upload Walkthrough

This guide details step-by-step procedures for different aspects of the GDC Data Submission process and how they relate to the GDC Data Model and structure. The first sections of this guide break down the submission process and associate each step with the Data Model. Additional sections are detailed below for strategies on expediting data submission, using features of the GDC Data Submission Portal, and best practices used by the GDC.

## GDC Data Model Basics

Pictured below is the submittable subset of the GDC Data Model: a roadmap for GDC data submission. Each entity type is represented with an oval in the graphic. All submitted entities require a connection to another entity type, based on the GDC Data Model, and a `submitter_id` as an identifier. This walkthrough will go through the submission of different entities. The completed (submitted) portion of the entity process will be highlighted in __blue__. 

[![GDC Data Model 1](images/GDC-Data-Model-None.png)](images/GDC-Data-Model-None.png "Click to see the full image.")

# Case and Clinical Data Submission

The `case` is the center of the GDC Data Model and usually describes a specific patient. Each `case` is connected to a `project`.  Different types of clinical data, such as `diagnoses` and `exposures`, are connected to the `case` to describe the case's attributes and medical information.   

## Case Submission

[![GDC Data Model 2](images/GDC-Data-Model-Case.png)](images/GDC-Data-Model-Case.png "Click to see the full image.")

The main entity of the GDC Data Model is the `case`, each of which must be registered beforehand with [dbGaP](https://www.ncbi.nlm.nih.gov/sra/docs/submitdbgap) under a unique `submitter_id`. The first step to submitting a `case` is to consult the [Data Dictionary](https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#data-dictionary-viewer), which details the fields that are associated with a `case`, the fields that are required to submit a `case`, and the values that can populate each field. Dictionary entries are available for all entities in the GDC Data Model.

[![Dictionary Case](images/Dictionary_Case.png)](images/Dictionary_Case.png "Click to see the full image.")

Submitting a [__Case__](https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=case) entity requires:

* __`submitter_id`:__ A unique key to identify the `case`
* __`projects.code`:__ A link to the `project`

The submitter ID is different from the universally unique identifier (UUID), which is based on the [UUID Version 4 Naming Convention](https://en.wikipedia.org/wiki/Universally_unique_identifier#Version_4_.28random.29). The UUID can be accessed under the `<entity_type>_id` field for each entity. For example, the `case` UUID can be accessed under the `case_id` field. The UUID is either assigned to each entity automatically or can be submitted by the user. Submitter-generated UUIDs cannot be uploaded in `submittable_data_file` entity types. See the [Data Model Users Guide](https://docs.gdc.cancer.gov/Data/Data_Model/GDC_Data_Model/#gdc-identifiers) for more details about GDC identifiers.

The `projects.code` field connects the `case` entity to the `project` entity.  The rest of the entity connections use the `submitter_id` field instead.  

The `case` entity can be added in JSON or TSV format. A template for any entity in either of these formats can be found in the Data Dictionary at the top of each page. Templates populated with `case` metadata in both formats are displayed below.  

```JSON
{
    "type": "case",
    "submitter_id": "PROJECT-INTERNAL-000055",
    "projects": {
        "code": "INTERNAL"
    }
}
```
```TSV
type  submitter_id  projects.code
case  PROJECT-INTERNAL-000055 INTERNAL   
```

>__Note:__ JSON and TSV formats handle links between entities (`case` and `project`) differently.  JSON includes the `code` field nested within `projects` while TSV appends `code` to `projects` with a period.  


## Uploading the Case Submission File

The file detailed above can be uploaded using the GDC Data Submission Portal and the GDC API as described below:

### Upload Using the GDC Data Submission Portal

An example of a `case` upload is detailed below. The [GDC Data Submission Portal](https://gdc.cancer.gov/submit-data/gdc-data-submission-portal) is equipped with a wizard window to facilitate the upload and validation of entities. 

#### 1. Upload Files

Choosing _'UPLOAD'_ from the project dashboard will open the Upload Data Wizard.

[![GDC Submission Wizard Upload Files](images/GDC_Submission_Wizard_Upload_2.png)](images/GDC_Submission_Wizard_Upload_2.png "Click to see the full image.")

Files containing one or more entities can be added either by clicking on `CHOOSE FILE(S)` or using drag and drop. Files can be removed from the Upload Data Wizard by clicking on the garbage can icon that is displayed next to the file after the file is selected for upload.

#### 2. Validate Entities

The __Validate Entities__ stage acts as a safeguard against submitting incorrectly formatted data to the GDC Data Submission Portal. During the validation stage, the GDC API will validate the content of uploaded entities against the Data Dictionary to detect potential errors. Invalid entities will not be processed and must be corrected by the user and re-uploaded before being accepted. A validation error report provided by the system can be used to isolate and correct errors.

When the first file is added, the wizard will move to the Validate section and the user can continue to add files. When all files have been added, choosing `VALIDATE` will run a test to check if the entities are valid for submission.

[![GDC Submission Wizard Validate Files](images/GDC_Submission_Portal_Validate.png)](images/GDC_Submission_Portal_Validate.png "Click to see the full image.")

#### 3. Commit or Discard Files
If the upload contains valid entities, a new transaction will appear in the latest transactions panel with the option to `COMMIT` or `DISCARD` the data. Entities contained in these files can be committed (applied) to the project or discarded using these two buttons.

If the upload contains invalid files, a transaction will appear with a FAILED status. Invalid files will need to be either corrected and re-uploaded or removed from the submission. If more than one file is uploaded and at least one is not valid, the validation step will fail for all files.  

[![Commit_Discard](images/GDC_Submission_CommitDiscard.png)](images/GDC_Submission_CommitDiscard.png "Click to see the full image.")


### Upload Using the GDC API

The API has a much broader range of functionality than the Data Wizard. Entities can be created, updated, and deleted through the API. See the [API Submission User Guide](https://docs.gdc.cancer.gov/API/Users_Guide/Submission/#creating-and-updating-entities) for a more detailed explanation and for the rest of the functionalities of the API. Generally, uploading an entity through the API can be performed using a command similar to the following:

```Shell
curl --header "X-Auth-Token: $token" --request POST --data @CASE.json https://api.gdc.cancer.gov/v0/submission/GDC/INTERNAL/_dry_run?async=true
```
CASE.json is detailed below.
```json
{
    "type": "case",
    "submitter_id": "PROJECT-INTERNAL-000055",
    "projects": {
        "code": "INTERNAL"
    }
}
```

In this example, the `_dry_run` marker is used to determine if the entities can be validated, but without committing any information. If a command passed through the `_dry_run` works, the command will work when it is changed to `commit`. For more information please go to [Dry Run Transactions](https://docs.gdc.cancer.gov/API/Users_Guide/Submission/#dry-run-transactions).

>__Note:__ Submission of TSV files is also supported by the GDC API.

Next, the file can either be committed (applied to the project) through the Data Submission Portal as before, or another API query can be performed that will commit the file to the project. The transaction number in the URL (467) is printed to the console during the first step of API submission and can also be retrieved from the [Transactions](Data_Submission_Process.md#transactions) tab in the Data Submission Portal.

```Shell
curl --header "X-Auth-Token: $token" --request POST https://api.gdc.cancer.gov/v0/submission/GDC/INTERNAL/transactions/467/commit?async=true
```

## Clinical Submission

Typically a submission project will include additional information about a `case` such as `demographic`, `diagnosis`, or `exposure` data.

### Clinical Data Requirements

For the GDC to release a project there is a minimum number of clinical properties that are required.  Minimal GDC requirements for each project includes age, gender, and diagnosis information.  Other [requirements](https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#?view=table-entity-list&anchor=clinical) may be added when the submitter is approved for submission to the GDC.

[![GDC Data Model Clinical](images/GDC-Data-Model-Clinical.png)](images/GDC-Data-Model-Clinical.png "Click to see the full image.")

### Submitting a Demographic Entity to a Case

The `demographic` entity contains information that characterizes the `case` entity.  

Submitting a [__Demographic__](https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=demographic) entity requires:

* __`submitter_id`:__ A unique key to identify the `demographic` entity.
* __`cases.submitter_id`:__ The unique key that was used for the `case` that links the `demographic` entity to the `case`.
* __`ethnicity`:__ An individual's self-described social and cultural grouping, specifically whether an individual describes themselves as Hispanic or Latino. The provided values are based on the categories defined by the U.S. Office of Management and Business and used by the U.S. Census Bureau.
* __`gender`:__ Text designations that identify gender. Gender is described as the assemblage of properties that distinguish people on the basis of their societal roles.
* __`race`:__ An arbitrary classification of a taxonomic group that is a division of a species. It usually arises as a consequence of geographical isolation within a species and is characterized by shared heredity, physical attributes and behavior, and in the case of humans, by common history, nationality, or geographic distribution. The provided values are based on the categories defined by the U.S. Office of Management and Business and used by the U.S. Census Bureau.
* __`year_of_birth`:__ Numeric value to represent the calendar year in which an individual was born.

```JSON
{
    "type": "demographic",
    "submitter_id": "PROJECT-INTERNAL-000055-DEMOGRAPHIC-1",
    "cases": {
        "submitter_id": "PROJECT-INTERNAL-000055"
    },
    "ethnicity": "not hispanic or latino",
    "gender": "male",
    "race": "asian",
    "year_of_birth": 1946
}
```
```TSV
type	cases.submitter_id	ethnicity	gender	race	year_of_birth
demographic	PROJECT-INTERNAL-000055	not hispanic or latino	male	asian	1946
```

### Submitting a Diagnosis Entity to a Case

Submitting a [__Diagnosis__](https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=diagnosis) entity requires:

* __`submitter_id`:__ A unique key to identify the `diagnosis` entity.
* __`cases.submitter_id`:__ The unique key that was used for the `case` that links the `diagnosis` entity to the `case`.
* __`age_at_diagnosis`:__ Age at the time of diagnosis expressed in number of days since birth.
* __`classification_of_tumor`:__ Text that describes the kind of disease present in the tumor specimen as related to a specific timepoint.
* __`days_to_last_follow_up`:__  Time interval from the date of last follow up to the date of initial pathologic diagnosis, represented as a calculated number of days.
* __`days_to_last_known_disease_status`:__ Time interval from the date of last follow up to the date of initial pathologic diagnosis, represented as a calculated number of days.
* __`days_to_recurrence`:__ Time interval from the date of new tumor event including progression, recurrence and new primary malignancies to the date of initial pathologic diagnosis, represented as a calculated number of days.
* __`last_known_disease_status`:__  The state or condition of an individual's neoplasm at a particular point in time.
* __`morphology`:__  The third edition of the International Classification of Diseases for Oncology, published in 2000 used principally in tumor and cancer registries for coding the site (topography) and the histology (morphology) of neoplasms. The study of the structure of the cells and their arrangement to constitute tissues and, finally, the association among these to form organs. In pathology, the microscopic process of identifying normal and abnormal morphologic characteristics in tissues, by employing various cytochemical and immunocytochemical stains. A system of numbered categories for representation of data.
* __`primary_diagnosis`:__  Text term for the structural pattern of cancer cells used to define a microscopic diagnosis.
* __`progression_or_recurrence`:__ Yes/No/Unknown indicator to identify whether a patient has had a new tumor event after initial treatment.
* __`site_of_resection_or_biopsy`:__ The third edition of the International Classification of Diseases for Oncology, published in 2000, used principally in tumor and cancer registries for coding the site (topography) and the histology (morphology) of neoplasms. The description of an anatomical region or of a body part. Named locations of, or within, the body. A system of numbered categories for representation of data.
* __`tissue_or_organ_of_origin`:__ Text term that describes the anatomic site of the tumor or disease.
* __`tumor_grade`:__ Numeric value to express the degree of abnormality of cancer cells, a measure of differentiation and aggressiveness.
* __`tumor_stage`:__ The extent of a cancer in the body. Staging is usually based on the size of the tumor, whether lymph nodes contain cancer, and whether the cancer has spread from the original site to other parts of the body. The accepted values for tumor_stage depend on the tumor site, type, and accepted staging system. These items should accompany the tumor_stage value as associated metadata.
* __`vital_status`:__ The survival state of the person registered on the protocol.

```JSON
{
    "type": "diagnosis",
    "submitter_id": "PROJECT-INTERNAL-000055-DIAGNOSIS-1",
    "cases": {
        "submitter_id": "GDC-INTERNAL-000099"
    },
    "age_at_diagnosis": 10256,
    "classification_of_tumor": "not reported",
    "days_to_last_follow_up": 34,
    "days_to_last_known_disease_status": 34,
    "days_to_recurrence": 45,
    "last_known_disease_status": "Tumor free",
    "morphology": "8260/3",
    "primary_diagnosis": "ACTH-producing tumor",
    "progression_or_recurrence": "no",
    "site_of_resection_or_biopsy": "Lung, NOS",
    "tissue_or_organ_of_origin": "Lung, NOS",
    "tumor_grade": "not reported",
    "tumor_stage": "stage i",
    "vital_status": "alive"
}
```
```TSV
type	submitter_id	cases.submitter_id	age_at_diagnosis	classification_of_tumor	days_to_last_follow_up	days_to_last_known_disease_status	days_to_recurrence	last_known_disease_status	morphology	primary_diagnosis	progression_or_recurrence	site_of_resection_or_biopsy	tissue_or_organ_of_origin	tumor_grade	tumor_stage	vital_status
diagnosis	PROJECT-INTERNAL-000055-DIAGNOSIS-1	GDC-INTERNAL-000099	10256	not reported	34	34	45	Tumor free	8260/3	ACTH-producing tumor	no	Lung, NOS	Lung, NOS	not reported	stage i	alive
```
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

### Submitting an Exposure Entity to a Case

Submitting an [__Exposure__](https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=exposure) entity does not require any information besides a link to the `case` and a `submitter_id`.  The following fields are optionally included:  

* __`alcohol_history`:__ A response to a question that asks whether the participant has consumed at least 12 drinks of any kind of alcoholic beverage in their lifetime.
* __`alcohol_intensity`:__ Category to describe the patient's current level of alcohol use as self-reported by the patient.
* __`bmi`:__ The body mass divided by the square of the body height expressed in units of kg/m^2.
* __`cigarettes_per_day`:__ The average number of cigarettes smoked per day (number).
* __`height`:__ The height of the individual in cm (number).
* __`weight`:__ The weight of the individual in kg (number).
* __`years_smoked`:__ Numeric value (or unknown) to represent the number of years a person has been smoking.

```JSON
{
    "type": "exposure",
    "submitter_id": "PROJECT-INTERNAL-000055-EXPOSURE-1",
    "cases": {
        "submitter_id": "PROJECT-INTERNAL-000055"
    },
    "alcohol_history": "yes",
    "bmi": 27.5,
    "cigarettes_per_day": 20,
    "height": 190,
    "weight": 100,
    "years_smoked": 5
}
```
```TSV
type	submitter_id	cases.submitter_id	alcohol_history	bmi	cigarettes_per_day	height	weight	years_smoked
exposure	PROJECT-INTERNAL-000055-EXPOSURE-1	PROJECT-INTERNAL-000055	yes	27.5	20	190	100	5
```

>__Note:__ Submitting a clinical entity uses the same conventions as submitting a `case` entity (detailed above).


# Biospecimen Submission

One of the main features of the GDC is the genomic data harmonization workflow. Genomic data is connected the case through biospecimen entities.  The `sample` entity describes a biological piece of matter that originated from a `case`.  Subsets of the `sample` such as `portions` and `analytes` can optionally be described.  The `aliquot` originates from a `sample` or `analyte` and describes the nucleic acid extract that was sequenced. The `read_group` entity describes the resulting set of reads from one sequencing lane.

## Sample Submission

[![GDC Data Model 3](images/GDC-Data-Model-Sample.png)](images/GDC-Data-Model-Sample.png "Click to see the full image.")

A `sample` submission has the same general structure as a `case` submission as it will require a unique key and a link to the `case`.  However, `sample` entities require one additional value:  `sample_type`. This peripheral data is required because it is necessary for the data to be interpreted. For example, an investigator using this data would need to know whether the `sample` came from tumor or normal tissue.  


[![Dictionary Sample](images/Dictionary_Sample.png)](images/Dictionary_Sample.png "Click to see the full image.")

Submitting a [__Sample__](https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=sample) entity requires:

* __`submitter_id`:__ A unique key to identify the `sample`.
* __`cases.submitter_id`:__ The unique key that was used for the `case` that links the `sample` to the `case`.
* __`sample_type`:__ Type of the `sample`. Named for its cellular source, molecular composition, and/or therapeutic treatment.

>__Note:__ The `case` must be "committed" to the project before a `sample` can be linked to it.  This also applies to all other links between entities.

```JSON
{
    "type": "sample",
    "cases": {
        "submitter_id": "PROJECT-INTERNAL-000055"
    },
    "sample_type": "Blood Derived Normal",
    "submitter_id": "Blood-00001SAMPLE_55"
}
```
```TSV
type	cases.submitter_id	submitter_id	sample_type
sample	PROJECT-INTERNAL-000055	Blood-00001SAMPLE_55	Blood Derived Normal  
```

## Portion, Analyte and Aliquot Submission

[![GDC Data Model 4](images/GDC-Data-Model-Aliquot.png)](images/GDC-Data-Model-Aliquot.png "Click to see the full image.")

Submitting a [__Portion__](https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=portion) entity requires:

* __`submitter_id`:__ A unique key to identify the `portion`.
* __`samples.submitter_id`:__ The unique key that was used for the `sample` that links the `portion` to the `sample`.

```JSON
{
    "type": "portion",
    "submitter_id": "Blood-portion-000055",
    "samples": {
        "submitter_id": "Blood-00001SAMPLE_55"
    }
}

```
```TSV
type	submitter_id	samples.submitter_id
portion	Blood-portion-000055	Blood-00001SAMPLE_55
```

Submitting an [__Analyte__](https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=analyte) entity requires:

* __`submitter_id`:__ A unique key to identify the `analyte`.
* __`portions.submitter_id`:__ The unique key that was used for the `portion` that links the `analyte` to the `portion`.
* __`analyte_type`:__ Protocol-specific molecular type of the specimen.

```JSON
{
    "type": "analyte",
    "portions": {
        "submitter_id": "Blood-portion-000055"
    },
    "analyte_type": "DNA",
    "submitter_id": "Blood-analyte-000055"
}

```
```TSV
type	portions.submitter_id	analyte_type	submitter_id
analyte	Blood-portion-000055	DNA	Blood-analyte-000055
```

Submitting an [__Aliquot__](https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=aliquot) entity requires:

* __`submitter_id`:__ A unique key to identify the `aliquot`.
* __`analytes.submitter_id`:__ The unique key that was used for the `analyte` that links the `aliquot` to the `analyte`.

```JSON
{
    "type": "aliquot",
    "submitter_id": "Blood-00021-aliquot55",
    "analytes": {
        "submitter_id": "Blood-analyte-000055"
    }
}

```
```TSV
type	submitter_id	analytes.submitter_id
aliquot	Blood-00021-aliquot55	Blood-analyte-000055
```

>__Note:__ `aliquot` entities can be directly linked to `sample` entities via the `samples.submitter_id`. The `portion` and `analyte` entities are not required for submission.

## Read Group Submission

[![GDC Data Model 5](images/GDC-Data-Model-RG.png)](images/GDC-Data-Model-RG.png "Click to see the full image.")

Information about sequencing reads is necessary for downstream analysis, thus the `read_group` entity requires more fields than the other Biospecimen entities (`sample`, `portion`, `analyte`, `aliquot`).

Submitting a [__Read Group__](https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=read_group) entity requires:

* __`submitter_id`:__ A unique key to identify the `read_group`.
* __`aliquots.submitter_id`:__ The unique key that was used for the `aliquot` that links the `read_group` to the `aliquot`.
* __`experiment_name`:__ Submitter-defined name for the experiment.
* __`is_paired_end`:__ Are the reads paired end? (Boolean value: `true` or `false`).
* __`library_name`:__ Name of the library.
* __`library_strategy`:__ Library strategy.
* __`platform`:__ Name of the platform used to obtain data.
* __`read_group_name`:__ The name of the `read_group`.
* __`read_length`:__ The length of the reads (integer).
* __`sequencing_center`:__ Name of the center that provided the sequence files.
* __`library_selection`:__ Library Selection Method.
* __`target_capture_kit`:__ Description that can uniquely identify a target capture kit. Suggested value is a combination of vendor, kit name, and kit version.

```JSON
{
    "type": "read_group",
    "submitter_id": "Blood-00001-aliquot_lane1_barcodeACGTAC_55",
    "experiment_name": "Resequencing",
    "is_paired_end": true,
    "library_name": "Solexa-34688",
    "library_strategy": "WXS",
    "platform": "Illumina",
    "read_group_name": "205DD.3-2",
    "read_length": 75,
    "sequencing_center": "BI",
    "library_selection": "Hybrid Selection",
    "target_capture_kit": "Custom MSK IMPACT Panel - 468 Genes",
    "aliquots":
        {
            "submitter_id": "Blood-00021-aliquot55"
        }    
}

```
```TSV
type	submitter_id	experiment_name	is_paired_end	library_name	library_strategy	platform	read_group_name	read_length	sequencing_center library_selection target_capture_kit	aliquots.submitter_id
read_group	Blood-00001-aliquot_lane1_barcodeACGTAC_55	Resequencing	true	Solexa-34688	WXS	Illumina	205DD.3-2	75	BI	Hybrid Selection Custom MSK IMPACT Panel - 468 Genes Blood-00021-aliquot55
```

>__Note:__ Submitting a biospecimen entity uses the same conventions as submitting a `case` entity (detailed above).

# Experiment Data Submission

Several types of experiment data can be uploaded to the GDC.  The `submitted_aligned_reads` and `submitted_unaligned_reads` files are associated with the `read_group` entity. While the array-based files such as the `submitted_tangent_copy_number` are associated with the `aliquot` entity.  Each of these file types are described in their respective entity submission and are uploaded separately using the [GDC API](https://docs.gdc.cancer.gov/API/Users_Guide/Getting_Started/) or the [GDC Data Transfer Tool](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool).  

[![GDC Data Model 6](images/GDC-Data-Model-Reads.png)](images/GDC-Data-Model-Reads.png "Click to see the full image.")

Before the experiment data file can be submitted, the GDC requires that the user provides information about the file as a `submittable_data_file` entity. This includes file-specific data needed to validate the file and assess which analyses should be performed. Sequencing data files can be submitted as `submitted_aligned_reads` or `submitted_unaligned_reads`.

Submitting a [__Submitted Aligned-Reads__](https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=submitted_aligned_reads) ([__Submitted Unaligned-Reads__](https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=submitted_unaligned_reads)) entity requires:

* __`submitter_id`:__ A unique key to identify the `submitted_aligned_reads`.
* __`read_groups.submitter_id`:__ The unique key that was used for the `read_group` that links the `submitted_aligned_reads` to the `read_group`.
* __`data_category`:__ Broad categorization of the contents of the data file.
* __`data_format`:__ Format of the data files.
* __`data_type`:__ Specific content type of the data file. (must be "Aligned Reads").
* __`experimental_strategy`:__ The sequencing strategy used to generate the data file.
* __`file_name`:__ The name (or part of a name) of a file (of any type).
* __`file_size`:__ The size of the data file (object) in bytes.
* __`md5sum`:__ The 128-bit hash value expressed as a 32 digit hexadecimal number used as a file's digital fingerprint.


```JSON
{
    "type": "submitted_aligned_reads",
    "submitter_id": "Blood-00001-aliquot_lane1_barcodeACGTAC_55.bam",
    "data_category": "Raw Sequencing Data",
    "data_format": "BAM",
    "data_type": "Aligned Reads",
    "experimental_strategy": "WGS",
    "file_name": "test.bam",
    "file_size": 38,
    "md5sum": "aa6e82d11ccd8452f813a15a6d84faf1",
    "read_groups": [
        {
            "submitter_id": "Primary_Tumor_RG_86-1"
        }
    ]
}
```
```TSV
type	submitter_id	data_category	data_format	data_type	experimental_strategy	file_name	file_size	md5sum	read_groups.submitter_id#1
submitted_aligned_reads	Blood-00001-aliquot_lane1_barcodeACGTAC_55.bam	Raw Sequencing Data	BAM	Aligned Reads	WGS	test.bam	38	aa6e82d11ccd8452f813a15a6d84faf1	Primary_Tumor_RG_86-1
```

>__Note:__ For details on submitting experiment data associated with more than one `read_group` entity, see the [Tips for Complex Submissions](#submitting-complex-data-model-relationships) section.    

## Uploading the Submittable Data File to the GDC

The submittable data file can be uploaded when it is registered with the GDC. An submittable data file is registered when its corresponding entity (e.g. `submitted_unaligned_reads`) is uploaded and committed. Uploading the file can be performed with either the [GDC Data Transfer Tool](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool) or the [GDC API](https://docs.gdc.cancer.gov/API/Users_Guide/Getting_Started/). Other types of data files such as clinical supplements, biospecimen supplements, and pathology reports are uploaded to the GDC in the same way. Supported data file formats are listed at the GDC [Submitted Data Types and File Formats](https://gdc.cancer.gov/about-data/data-types-and-file-formats/submitted-data-types-and-file-formats) website.

__GDC Data Transfer Tool:__ A file can be uploaded using its UUID (which can be retrieved from the GDC Submission Portal or API) once it is registered. 

[![UUID Location](images/GDC_Submission_UUID_location.png)](images/GDC_Submission_UUID_location.png "Click to see the full image.")

The following command can be used to upload the file:

```Shell
gdc-client upload --project-id PROJECT-INTERNAL --identifier a053fad1-adc9-4f2d-8632-923579128985 -t $token -f $path_to_file
```   

Additionally a manifest can be downloaded from the Submission Portal and passed to the Data Transfer Tool. This will allow for the upload of more than one `submittable_data_file`:

```Shell
gdc-client upload -m manifest.yml -t $token
```
__API Upload:__  A `submittable_data_file` can be uploaded through the API by using the `/submission/$PROGRAM/$PROJECT/files` endpoint.  The following command would be typically used to upload a file:  

```Shell
curl --request PUT --header "X-Auth-Token: $token" https://api.gdc.cancer.gov/v0/submission/PROJECT/INTERNAL/files/6d45f2a0-8161-42e3-97e6-e058ac18f3f3 -d $path_to_file

```

For more details on how to upload a `submittable_data_file` to a project see the [API Users Guide](https://docs.gdc.cancer.gov/API/Users_Guide/Submission/) and the [Data Transfer Tool Users Guide](https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Data_Download_and_Upload/).  

## Annotation Submission

The GDC Data Portal supports the use of annotations for any submitted entity or file.  An annotation entity may include comments about why particular patients or samples are not present or why they may exhibit critical differences from others.  Annotations include information that cannot be submitted to the GDC through other existing nodes or properties.

If a submitter would like to create an annotation please contact the GDC Support Team (support@nci-gdc.datacommons.io).

## Deleting Submitted Entities

The GDC Data Submission Portal allows users to delete submitted entities from the project when the project is in an "OPEN" state. Files cannot be deleted while in the "SUBMITTED" state.  This section applies to entities that have been committed to the project. Entities that have not been committed can be removed from the project by choosing the `DISCARD` button.  Entities can also be deleted using the API. See the [API Submission Documentation](https://docs.gdc.cancer.gov/API/Users_Guide/Submission/#deleting-entities) for specific instructions.

>__NOTE:__  Entities associated with files uploaded to the GDC object store cannot be deleted until the associated file has been deleted. Users must utilize the [GDC Data Transfer Tool](https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Data_Download_and_Upload/#deleting-previously-uploaded-data) to delete these files first.

### Simple Deletion

If an entity was uploaded and has no related entities, it can be deleted from the [Browse](Data_Submission_Process.md#browse-data) tab. Once the entity to be deleted is selected, choose the `DELETE` button in the right panel under "ACTIONS".


[![GDC Delete Unassociated Case](images/GDC-Delete-Case-Unassociated.png)](images/GDC-Delete-Case-Unassociated.png "Click to see the full image.")


A message will then appear asking if you are sure about deleting the entity.  Choosing the `YES, DELETE` button will remove the entity from the project, whereas choosing the `NO, CANCEL` button will return the user to the previous screen.  


[![GDC Yes or No](images/GDC-Delete-Sure.png)](images/GDC-Delete-Sure.png "Click to see the full image.")


### Deletion with Dependents

If an entity has related entities, such as a `case` with multiple `samples` and `aliquots`, deletion takes one extra step.  


[![GDC Delete Associated Case](images/GDC-Delete-Case-Associated.png)](images/GDC-Delete-Case-Associated.png "Click to see the full image.")


Follow the [Simple Deletion](Data_Submission_Walkthrough.md#simple-deletion) method until the end. This action will appear in the [Transactions](Data_Submission_Process.md#transactions) tab as "Delete" with a "FAILED" state.  


[![GDC Delete Failed](images/GDC-Failed-Transaction.png)](images/GDC-Failed-Transaction.png "Click to see the full image.")


Choose the failed transaction and the right panel will show the list of entities related to the entity that was going to be deleted.  


[![GDC Error Related](images/GDC-Error-Related.png)](images/GDC-Error-Related.png "Click to see the full image.")


Selecting the `DELETE ALL` button at the bottom of the list will delete all of the related entities, their descendants, and the original entity.


### Submitted Data File Deletion

The [`submittable_data_files`](https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#?view=table-entity-list&anchor=submittable_data_file) that were uploaded erroneously are deleted separately from their associated entity using the GDC Data Transfer Tool. See the section on [Deleting Data Files](https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Data_Download_and_Upload/#deleting-previously-uploaded-data) in the Data Transfer Tool users guide for specific instructions.  

## Updating Uploaded Entities

Before harmonization occurs, entities can be modified to update, add, or delete information. These methods are outlined below.

### Updating or Adding Fields

Updated or additional fields can be applied to entities by re-uploading them through the GDC Data Submission portal or API. See below for an example of a case upload with a `primary_site` field being added and a `disease_type` field being updated.

```Before
{
"type":"case",
"submitter_id":"GDC-INTERNAL-000043",
"projects":{
  "code":"INTERNAL"
},
"disease_type": "Neuroblastoma"
}
```
```After
{
"type":"case",
"submitter_id":"GDC-INTERNAL-000043",
"projects":{
  "code":"INTERNAL"
},
"disease_type": "Germ Cell Neoplasms",
"primary_site": "Pancreas"
}
```
__Guidelines:__

* The newly uploaded entity must contain the `submitter_id` of the existing entity so that the system updates the correct one.
* All newly updated entities will be validated by the GDC Dictionary.  All required fields must be present in the newly updated entity.
* Fields that are not required do not need to be re-uploaded and will remain unchanged in the entity unless they are updated.

### Deleting Optional Fields

It may be necessary to delete fields from uploaded entities. This can be performed through the API and can only be applied to optional fields. It also requires the UUID of the entity, which can be retrieved from the submission portal or using a GraphQL query.

In the example below, the `primary_site` and `disease_type` fields are removed from a `case` entity:

```Shell
curl --header "X-Auth-Token: $token_string" --request DELETE  --header "Content-Type: application/json" "https://api.gdc.cancer.gov/v0/submission/EXAMPLE/PROJECT/entities/7aab7578-34ff-5651-89bb-57aefdc4c4f8?fields=primary_site,disease_type"
```

```Before
{
"type":"case",
"submitter_id":"GDC-INTERNAL-000043",
"projects":{
  "code":"INTERNAL"
},
"disease_type": "Germ Cell Neoplasms",
"primary_site": "Pancreas"
}
```
```After
{
"type":"case",
"submitter_id":"GDC-INTERNAL-000043",
"projects":{
  "code":"INTERNAL"
}
}
```

### Versioning
Changes to entities will create versions. For more information on this, please go to [Uploading New Versions of Data Files](https://docs.gdc.cancer.gov/API/Users_Guide/Submission/#uploading-new-versions-of-data-files).

## Strategies for Submitting in Bulk

Each submission in the previous sections was broken down by component to demonstrate the GDC Data Model structure. However, the submission of multiple entities at once is supported and encouraged. Here two strategies for submitting data in an efficient manner are discussed.   

### Registering a BAM File: One Step

Registering a BAM file (or any other type) can be performed in one step by including all of the entities, from `case` to `submitted_aligned_reads`, in one file.  See the example below:

```JSON
[{
    "type": "case",
    "submitter_id": "PROJECT-INTERNAL-000055",
    "projects": {
        "code": "INTERNAL"
    }
},
{
    "type": "sample",
    "cases": {
        "submitter_id": "PROJECT-INTERNAL-000055"
    },
    "sample_type": "Blood Derived Normal",
    "submitter_id": "Blood-00001_55"
},
{
    "type": "portion",
    "submitter_id": "Blood-portion-000055",
    "samples": {
        "submitter_id": "Blood-00001_55"
    }
},
{
    "type": "analyte",
    "portions": {
        "submitter_id": "Blood-portion-000055"
    },
    "analyte_type": "DNA",
    "submitter_id": "Blood-analyte-000055"
},
{
    "type": "aliquot",
    "submitter_id": "Blood-00021-aliquot55",
    "analytes": {
        "submitter_id": "Blood-analyte-000055"
    }
},
{
    "type": "read_group",
    "submitter_id": "Blood-00001-aliquot_lane1_barcodeACGTAC_55",
    "experiment_name": "Resequencing",
    "is_paired_end": true,
    "library_name": "Solexa-34688",
    "library_strategy": "WXS",
    "platform": "Illumina",
    "read_group_name": "205DD.3-2",
    "read_length": 75,
    "sequencing_center": "BI",
    "aliquots":
        {
            "submitter_id": "Blood-00021-aliquot55"
        }    
},
{
    "type": "submitted_aligned_reads",
    "submitter_id": "Blood-00001-aliquot_lane1_barcodeACGTAC_55.bam",
    "data_category": "Raw Sequencing Data",
    "data_format": "BAM",
    "data_type": "Aligned Reads",
    "experimental_strategy": "WGS",
    "file_name": "test.bam",
    "file_size": 38,
    "md5sum": "aa6e82d11ccd8452f813a15a6d84faf1",
    "read_groups": [
        {
            "submitter_id": "Blood-00001-aliquot_lane1_barcodeACGTAC_55"
        }
    ]
}]
```

All of the entities are placed into a JSON list object:

`[{"type": "case","submitter_id": "PROJECT-INTERNAL-000055","projects": {"code": "INTERNAL"}}}, entity-2, entity-3]`

The entities need not be in any particular order as they are validated together.

>__Note:__ Tab-delimited format is not recommended for 'one-step' submissions due to an inability of the format to accommodate multiple 'types' in one row.  

### Submitting Numerous Cases

The GDC understands that submitters will have projects that comprise more entities than would be reasonable to individually parse into JSON formatted files. Additionally, many investigators store large amounts of data in a tab-delimited format (TSV).  For instances like this, we recommend parsing all entities of the same type into separate TSVs and submitting them on a type-basis.  

For example, a user may want to submit 100 Cases associated with 100 `samples`, 100 `portions`, 100 `analytes`, 100 `aliquots`, and 100 `read_groups`. Constructing and submitting 100 JSON files would be tedious and difficult to organize. The solution is submitting one `case` TSV containing the 100 `cases`, one `sample` TSV containing the 100 `samples`, so on and so forth. Doing this would only require six TSVs and these files can be formatted in programs such as Microsoft Excel or Google Spreadsheets.  

See the following example TSV files:

* [Cases.tsv](Cases.tsv)
* [Samples.tsv](Samples.tsv)
* [Portions.tsv](Portions.tsv)
* [Analytes.tsv](Analytes.tsv)
* [Aliquots.tsv](Aliquots.tsv)
* [Read-Groups.tsv](Readgroups.tsv)

### Download Previously Uploaded Metadata Files

The [transaction](Data_Submission_Process.md#transactions) page lists all previous transactions in the project. The user can download metadata files uploaded to the GDC workspace in the details section of the screen by selecting one transaction and scrolling to the "DOCUMENTS" section.


[![Transaction Original Files](images/GDC_Submission_Transactions_Original_Files_2.png)](images/GDC_Submission_Transactions_Original_Files_2.png "Click to see the full image.")

### Download Previously Uploaded Data Files

The only supported method to download data files previously uploaded to the GDC Submission Portal that have not been release yet is to use the API or the [Data Transfer Tool](https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Getting_Started/). To retrieve data previous upload to the submission portal you will need to retrieve the data file's UUID.  The UUIDs for submitted data files are located in the submission portal under the file's Summary section as well as the manifest file located on the file's Summary page. 

[![Submission Portal Summary View](images/gdc-submission__image2_submission_UUID.png)](images/gdc-submission__image2_submission_UUID.png "Click to see the full image.")  

Once the UUID(s) have been retrieved, the download process is the same as it is for downloading data files at the [GDC Portal using UUIDs](https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Data_Download_and_Upload/#downloading-data-using-gdc-file-uuids).

 >__Note:__ When submittable data files are uploaded through the Data Transfer Tool they are not displayed as transactions.


# Tips for Complex Submissions

Because of the data types and relationships included in the GDC, data submission can become a complex procedure. The purpose of this section is to present guidelines that will aid in the incorporation and harmonization of submitters' data. Please contact the GDC Help Desk at __support@nci-gdc.datacommons.io __ if you have any questions or concerns regarding a submission project.


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

Note that we will only run variant calling for a particular tumor aliquot per experimental strategy once.  You must make sure that the appropriate normal control is uploaded to the GDC when Requesting Submission.  Uploading a different normal sample later will not result in reanalysis by the GDC.