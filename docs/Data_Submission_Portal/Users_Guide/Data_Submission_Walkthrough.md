# Data Upload Walkthrough

This guide details step-by-step procedures for different aspects of the GDC Data Submission process and how they relate to the GDC Data Model and structure. The first sections of this guide break down the submission process and associate each step with the Data Model. Additional sections are detailed below for strategies on expediting data submission, using features of the GDC Data Submission Portal, and best practices used by the GDC.

## GDC Data Model Basics

Pictured below is the submittable subset of the GDC Data Model: a roadmap for GDC data submission. Each oval node in the graphic represents an entity: a logical unit of data related to a specific clinical, biospecimen, or file facet in the GDC. An entity includes a set of fields, the associated values, and information about its related node associations. All submitted entities require a connection to another entity type, based on the GDC Data Model, and a `submitter_id` as an identifier. This walkthrough will go through the submission of different entities. The completed (submitted) portion of the entity process will be highlighted in __blue__. 

[![GDC Data Model 1](images/GDC-Data-Model-None.png)](images/GDC-Data-Model-None.png "Click to see the full image.")

# Case Submission

The `case` is the center of the GDC Data Model and usually describes a specific patient. Each `case` is connected to a `project`.  Different types of clinical data, such as `diagnoses` and `exposures`, are connected to the `case` to describe the case's attributes and medical information.   

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

# Clinical Data Submission

Typically, a submission project will include additional information about a `case` such as `demographic`, `diagnosis`, or `exposure` data.

## Clinical Data Requirements

For the GDC to release a project there is a minimum number of clinical properties that are required.  Minimal GDC requirements for each project includes age, gender, and diagnosis information.  Other [requirements](https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#?view=table-entity-list&anchor=clinical) may be added when the submitter is approved for submission to the GDC.

[![GDC Data Model Clinical](images/GDC-Data-Model-Clinical.png)](images/GDC-Data-Model-Clinical.png "Click to see the full image.")

## Submitting a Demographic Entity to a Case

The `demographic` entity contains information that characterizes the `case` entity.  

Submitting a [__Demographic__](https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=demographic) entity requires:

* __`submitter_id`:__ A unique key to identify the `demographic` entity.
* __`cases.submitter_id`:__ The unique key that was used for the `case` that links the `demographic` entity to the `case`.
* __`ethnicity`:__ An individual's self-described social and cultural grouping, specifically whether an individual describes themselves as Hispanic or Latino. The provided values are based on the categories defined by the U.S. Office of Management and Business and used by the U.S. Census Bureau.
* __`gender`:__ Text designations that identify gender. Gender is described as the assemblage of properties that distinguish people on the basis of their societal roles.
* __`race`:__ An arbitrary classification of a taxonomic group that is a division of a species. It usually arises as a consequence of geographical isolation within a species and is characterized by shared heredity, physical attributes and behavior, and in the case of humans, by common history, nationality, or geographic distribution. The provided values are based on the categories defined by the U.S. Office of Management and Business and used by the U.S. Census Bureau.

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
}
```
```TSV
type	cases.submitter_id	ethnicity	gender	race
demographic	PROJECT-INTERNAL-000055	not hispanic or latino	male	asian
```

## Submitting a Diagnosis Entity to a Case

Submitting a [__Diagnosis__](https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=diagnosis) entity requires:

* __`submitter_id`:__ A unique key to identify the `diagnosis` entity.
* __`cases.submitter_id`:__ The unique key that was used for the `case` that links the `diagnosis` entity to the `case`.
* __`age_at_diagnosis`:__ Age at the time of diagnosis expressed in number of days since birth.
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
type	submitter_id	cases.submitter_id	age_at_diagnosis    days_to_last_follow_up	days_to_last_known_disease_status	days_to_recurrence	last_known_disease_status	morphology	primary_diagnosis	progression_or_recurrence	site_of_resection_or_biopsy	tissue_or_organ_of_origin	tumor_grade	tumor_stage	vital_status
diagnosis	PROJECT-INTERNAL-000055-DIAGNOSIS-1	GDC-INTERNAL-000099	10256	34	34	45	Tumor free	8260/3	ACTH-producing tumor	no	Lung, NOS	Lung, NOS	not reported	stage i	alive
```

### Submitting an Exposure Entity to a Case

Submitting an [__Exposure__](https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=exposure) entity does not require any information besides a link to the `case` and a `submitter_id`.  The following fields are optionally included:  

* __`alcohol_history`:__ A response to a question that asks whether the participant has consumed at least 12 drinks of any kind of alcoholic beverage in their lifetime.
* __`alcohol_intensity`:__ Category to describe the patient's current level of alcohol use as self-reported by the patient.
* __`alcohol_days_per_week`:__ Numeric value used to describe the average number of days each week that a person consumes an alchoolic beverage.
* __`years_smoked`:__ Numeric value (or unknown) to represent the number of years a person has been smoking.
* __`tobacco_smoking_onset_year`:__ The year in which the participant began smoking.
* __`tobacco_smoking_quit_year`:__ The year in which the participant quit smoking. 

```JSON
{
    "type": "exposure",
    "submitter_id": "PROJECT-INTERNAL-000055-EXPOSURE-1",
    "cases": {
        "submitter_id": "PROJECT-INTERNAL-000055"
    },
    "alcohol_history": "yes",
    "alcohol_intensity": "Drinker",
    "alcohol_days_per_week": 2,
    "years_smoked": 5,
    "tobacco_smoking_onset_year": 2007,
    "tobacco_smoking_quit_year": 2012
}
```
```TSV
type	submitter_id	cases.submitter_id	alcohol_history alcohol_intensity   alcohol_days_per_week years_smoked  tobacco_smoking_onset_year  tobacco_smoking_quit_year
exposure	PROJECT-INTERNAL-000055-EXPOSURE-1	PROJECT-INTERNAL-000055	yes Drinker 2 5 2007    2012
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
* __`tissue_type`:__ Text term that represents a description of the kind of tissue collected with respect to disease status or proximity to tumor tissue.

>__Note:__ The `case` must be "committed" to the project before a `sample` can be linked to it.  This also applies to all other links between entities.

```JSON
{
    "type": "sample",
    "cases": {
        "submitter_id": "PROJECT-INTERNAL-000055"
    },
    "sample_type": "Blood Derived Normal",
    "submitter_id": "Blood-00001SAMPLE_55"
    "tissue_type": "Normal"
}
```
```TSV
type	cases.submitter_id	submitter_id	sample_type tissue_type
sample	PROJECT-INTERNAL-000055	Blood-00001SAMPLE_55	Blood Derived Normal    Normal
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
* __`analyte_type`:__ Text term that represents the kind of molecular specimen analyte.

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
type	submitter_id	experiment_name	is_paired_end	library_name    library_selection	library_strategy	platform	read_group_name	read_length	sequencing_center   target_capture_kit	aliquots.submitter_id
read_group	Blood-00001-aliquot_lane1_barcodeACGTAC_55	Resequencing	true	Solexa-34688    Hybrid Selection	WXS	Illumina	205DD.3-2	75	BI	Custom MSK IMPACT Panel - 468 Genes Blood-00021-aliquot55
```

>__Note:__ Submitting a biospecimen entity uses the same conventions as submitting a `case` entity (detailed above).

# Experiment Data Submission

Several types of experiment data can be uploaded to the GDC.  The `submitted_aligned_reads` and `submitted_unaligned_reads` files are associated with the `read_group` entity, while the array-based files such as the `submitted_tangent_copy_number` are associated with the `aliquot` entity.  Each of these file types are described in their respective entity submission and are uploaded separately using the [GDC API](https://docs.gdc.cancer.gov/API/Users_Guide/Getting_Started/) or the [GDC Data Transfer Tool](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool).  

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

The submittable data file can be uploaded when it is registered with the GDC. An submittable data file is registered when its corresponding entity (e.g. `submitted_unaligned_reads`) is uploaded and committed. It is important to not that the Harmonization process does not occur on these submitted files until the user clicks the [`Request Submission`](Data_Submission_Process.md#release) button. Uploading the file can be performed with either the [GDC Data Transfer Tool](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool) or the [GDC API](https://docs.gdc.cancer.gov/API/Users_Guide/Getting_Started/). Other types of data files such as clinical supplements, biospecimen supplements, and pathology reports are uploaded to the GDC in the same way. Supported data file formats are listed at the GDC [Submitted Data Types and File Formats](https://gdc.cancer.gov/about-data/data-types-and-file-formats/submitted-data-types-and-file-formats) website.

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

If a submitter would like to create an annotation, please contact the GDC Support Team (support@nci-gdc.datacommons.io).

## Deleting Submitted Entities

The GDC Data Submission Portal allows users to delete submitted entities from the project when the project is in an "OPEN" state. Files cannot be deleted while in the "SUBMITTED" state.  This section applies to entities that have been committed to the project. Entities that have not been committed can be removed from the project by choosing the `DISCARD` button.  Entities can also be deleted using the API. See the [API Submission Documentation](https://docs.gdc.cancer.gov/API/Users_Guide/Submission/#deleting-entities) for specific instructions.

>__NOTE:__  Entities associated with files uploaded to the GDC object store cannot be deleted until the associated file has been deleted. Users must utilize the [GDC Data Transfer Tool](https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Data_Download_and_Upload/#deleting-previously-uploaded-data) to delete these files first.

### Simple Deletion

If an entity was uploaded and has no related entities, it can be deleted from the [Browse](Data_Submission_Process.md#browse) tab. Once the entity to be deleted is selected, choose the `DELETE` button in the right panel under "ACTIONS".


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
"disease_type": "Myomatous Neoplasms"
}
```
```After
{
"type":"case",
"submitter_id":"GDC-INTERNAL-000043",
"projects":{
  "code":"INTERNAL"
},
"disease_type": "Myxomatous Neoplasms",
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
    "library_selection":"Hybrid Selection",
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
