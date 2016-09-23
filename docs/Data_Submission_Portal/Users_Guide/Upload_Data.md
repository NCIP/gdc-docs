# Upload and Validate Data

## Overview

The GDC Data Submission process is detailed on the [GDC Website]( https://gdc.nci.nih.gov/submit-data/data-submission-processes-and-tools).

This chapter will focus on the upload and validation of data to the project workspace.

[![GDC Data Submission Portal Workflow Upload](images/GDC_Submission_Portal_Workflow_Upload.png)](images/GDC_Submission_Portal_Workflow_Upload.png "Click to see the full image.")


## Introduction to Files

### File Format

The GDC Data Submission Portal supports upload for the following file formats:

* JSON
* TSV

During the upload and validation process, files are converted by the GDC API into entities and inserted into the database, maintaining a file-agnostic back end.

The GDC Data Submission Portal offers the ability to download files in different formats. The system converts database entities back to the requested file format.

### File Type
The GDC Data Submission Portal supports the following types of files for upload to the GDC:

* __Clinical__: A case's clinical data.
* __Biospecimen__: Metadata describing a tissue specimen collected from a case and other material derived from samples for analysis.
* __Annotations__: Observations associated with any entity, which can be useful for interpreting the data.

More details about the submission process, data files and file formats can be found on the GDC website at [Data Submission Processes and Tools](https://gdc.nci.nih.gov/submit-data/data-submission-processes-and-tools) and [Data Types and File Formats](https://gdc.nci.nih.gov/submit-data/gdc-data-types-and-file-formats).

Data supported by the GDC is described in the [GDC Data Dictionary](../../Data_Dictionary/viewer.md).

The table below is an example of files used to upload a read group to the GDC. The next section will describe how to perform these actions.

|# | File name | Description | GDC Tool |
| --- | --- | --- | --- |
|1a | read_group.tsv |Read Group Metadata|Upload this file to the Submission Portal. It describes the experiment metadata.|
|1b | submitted\_unaligned_reads.tsv|Submitted File Metadata|Upload this file to the Submission Portal. It describes the file metadata.|
|2 | manifest.yml|Manifest|Download this file from the Submission Portal. This manifest is used by the GDC Data Transfer Tool for the actual file upload.|
|3 | ExperimentFile.fastq|FASTQ File|Upload this file to the GDC Data Transfer Tool using the manifest (2).|



## Step 1: Prepare Files

The [GDC Data Dictionary](../../Data_Dictionary/viewer.md) describes the types of entities that can be uploaded to the GDC.

The user can go to the GDC Data Dictionary to __download the template files__ for a file. The templates can be populated with data by the user and should result in a valid upload (if validation rules detailed in the Data Dictionary are met). See each individual Data Dictionary entry for required fields and acceptable values for each field.


### Focus on Links

In order to identify the relationship between two entities, the user should include in the file of the child entity a reference to the parent submitter ID (called link).

For example, a Demographic entity describes a Case entity. The user is required to define __cases.submitter_id__ in the Demographic file.


### Examples

#### Demographic


An example of a __demographic file__ that can be uploaded in TSV format is detailed below.

The structure of the file is as follows:

* Type = demographic
* Unique Keys = project\_id, submitter_id
* Links = cases.submitter_id
* Properties = ethnicity, gender, etc.

```tsv
type	project_id	submitter_id	cases.submitter_id	ethnicity	gender	race	year_of_birth	year_of_death
demographic	TCGA-DEV3	TCGA-DEV-3-CASE-000-D1	TCGA-DEV-3-CASE-000	hispanic or latino	male	white	1950	0
demographic	TCGA-DEV3	TCGA-DEV-3-CASE-001-D1	TCGA-DEV-3-CASE-001	not reported	female	white	1956	0
```


#### Read Group and File

An example of a __Read Group__ upload, which requires two TSV files to describe metadata, is detailed below.

File 1: read_group.tsv

```tsv
type	project_id	submitter_id	aliquots.submitter_id	experiment_name	is_paired_end	library_name	library_strategy	platform	read_group_name	read_length	sequencing_center	RIN	adapter_name	adapter_sequence	base_caller_name	base_caller_version	fastq_name	flow_cell_barcode	includes_spike_ins	instrument_model	library_preparation_kit_catalog_number	library_preparation_kit_name	library_preparation_kit_vendor	library_preparation_kit_version	library_selection	library_strand	sequencing_date	size_selection_range	spike_ins_concentration	spike_ins_fasta	target_capture_kit_catalog_number	target_capture_kit_name	target_capture_kit_target_region	target_capture_kit_vendor	target_capture_kit_version	to_trim_adapter_sequence
read_group	TCGA-DEV3	read_group_ID1	TCGA-DEV-3-CASE-000-S1-AL1	Text for Experiment	TRUE	lib_1	WXS	Illumina	35	101	test								FALSE																	FALSE
```

File 2: submitted\_unaligned_reads.tsv

```tsv
type	project_id	submitter_id	read_groups.submitter_id	file_name	file_size	md5sum	data_category	data_type	data_format	experimental_strategy
submitted_unaligned_reads	TCGA-DEV3	fileID1_CASE-000-AL1	read_group_ID1	fileID88_CASE-000.fastq	61004	311253B0CA93B396A41C0A88F01557AE	Sequencing Data	Unaligned Reads	FASTQ	WGS
```


Once the files have been prepared (in TSV or JSON format), these files can be uploaded with the Upload Data Wizard.

**Note:** Before you can upload clinical, biospecimen or experiment data, associated cases must be registered in the GDC. If the cases are not displayed in your project dashboard, download the case template from the [GDC Data Dictionary](../../Data_Dictionary/viewer.md). Complete it with the Case Submitter IDs and upload the Cases with the Upload Data Wizard.

## Step 2: Upload Data Wizard

The GDC Data Submission Portal is equipped with a wizard window to guide you through the upload and validation of data. There are three stages:

* __Upload Files__: Upload a file into the user's browser, at this point nothing is submitted to the project workspace.
* __Validate Files__: Send a file to the GDC backend to validate its content (see below).
* __Confirm Upload__: Submit a validated file to the project workspace and produce a report.

The _'File Validation'_ stage acts as a safeguard against submitting incorrect files to the GDC Data Submission Portal. During the validation stage, the GDC API will validate the content of submitted files against the Data Dictionary to detect potential errors. Invalid files will be flagged and denied upload to the GDC until corrections are made by the user. A validation error report provided by the system can be used to isolate and correct errors for resubmission.

### Upload Files

From the project dashboard (panel 1), clicking on _'UPLOAD'_ will open the submission wizard.

[![GDC Submission Wizard Upload Files](images/GDC_Submission_Wizard_Upload_2.png)](images/GDC_Submission_Wizard_Upload_2.png "Click to see the full image.")

Files can be added either by clicking on _'CHOOSE FILE(S)'_ or by using drag and drop.

### Validate Files

As soon as the first file is added, the wizard will move to the _'VALIDATE'_ section and the user can continue to add files.

[![GDC Submission Wizard Validate Files](images/GDC_Submission_Wizard_Validate.png)](images/GDC_Submission_Wizard_Validate.png "Click to see the full image.")

Once all files have been added, clicking on _'VALIDATE'_ will check if the files are valid for submission.

[![Confirm Submission](images/GDC_Submission_Wizard_Confirm.png)](images/GDC_Submission_Wizard_Confirm.png "Click to see the full image.")

If the upload contains invalid files, the user will not be able to submit the data and those files will need to be either corrected and re-uploaded or removed from the submission.

[![Invalid Files in a Submission](images/GDC_Submission_Wizard_Invalid_Files.png)](images/GDC_Submission_Wizard_Invalid_Files.png "Click to see the full image.")

Files can be removed from the submission by clicking on the _'garbage can'_ icon related to the file.

Click on _'View Report'_ to see the error report which provides details about why files failed validation. This report is detailed in the [Data Submission Reports](Reports.md) section of the documentation.



### Confirm Upload

Once all files are valid for upload, clicking on _'Confirm Upload'_ will upload your files to the GDC project workspace.

**Note:** Uploaded files are not released immediately. The user must submit the data to the GDC for harmonization and release the project to make data available on the GDC Data Portal.

[![Successful Submission](images/GDC_Submission_wizard_successful_submission.png)](images/GDC_Submission_wizard_successful_submission.png "Click to see the full image.")

You can then click on _'CLICK HERE TO VIEW THE TRANSACTION'_ to be redirected to the transaction list page.

## Step 3: GDC Data Transfer Tool

The GDC Data Transfer Tool is used to upload the actual file.

Once the user has uploaded metadata through the Upload Data Wizard (e.g., read_group and submitted file metadata),  the manifest will be available for download from the transaction report.

**Note:** You can also download the manifest from the "Submittable Data Files" section of the Browse menu.

[![Transaction Manifest](images/GDC_Submission_Transactions_Get_Manifest_2.png)](images/GDC_Submission_Transactions_Get_Manifest_2.png "Click to see the full image.")

Users can use this manifest to upload their actual files with the GDC Data Transfer Tool. Please refer to the [GDC Data Transfer Tool's User Guide](../Data_Transfer_Tool/Users_Guide/Getting_Started.md) for more information.

## Download Previously Uploaded Files

The [transaction](Transactions.md) page lists all previous transactions in the project. The user can download files uploaded to the GDC workspace in the details section of the screen by selecting a particular transaction and scrolling to the "DOCUMENTS" section.

[![Transaction Original Files](images/GDC_Submission_Transactions_Original_Files_2.png)](images/GDC_Submission_Transactions_Original_Files_2.png "Click to see the full image.")
