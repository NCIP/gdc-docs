# Upload and Validate Data

## Overview

The GDC Data Submission process is detailed on the [GDC Website]( https://gdc.nci.nih.gov/submit-data/data-submission-processes-and-tools).

This chapter will focus on the upload and validation of data to the GDC project workspace.

## Introduction to the Files

### File Format

The GDC Data Submission Portal supports the following file formats for submission:

* JSON
* TSV

During the upload and validation process, files are converted by the GDC API into entities and inserted into the database, maintaining a file-agnostic back-end.

The GDC Data Submission Portal offers the ability to download files in different formats. To do so the system converts database entities back to the requested file format.

### File Type
The GDC Data Submission Portal supports four types of files for upload to the GDC:

* __Clinical__: A case’s clinical data.
* __Biospecimen__: Metadata describing a tissue specimen collected from a case and other material derived from samples for analysis.
* __Experiment Data__: GDC's unit of submission (see below).
* __Annotations__: Observations associated with any entity, which can be useful for interpreting the data.

More details about the submission process, data files and file formats can be found on the GDC website at [Data Submission Processes and Tools](https://gdc.nci.nih.gov/submit-data/data-submission-processes-and-tools) and [Data Types and File Formats](https://gdc.nci.nih.gov/submit-data/gdc-data-types-and-file-formats).

#### Focus on Experiment Data

The GDC has developed a submission unit called a “data bundle”, which is a set of files with associated metadata.

Data bundle types are defined by what files are expected and what introspection is done for validation or linking to other GDC entities. Each data bundle will be validated via a bundle type and project specific schema, including a JSON data dictionary, relationship check and molecular data quality check.

The experiment data supported by the GDC is described in the [GDC Dictionary](../../Dictionary/viewer.md).

The table below is an example of files used to upload a read group to the GDC. The next section will describe how to perform these actions.

| File | Example | Usage |
| --- | --- | --- |
| Read Group Metadata|read_group.tsv|Upload this file to the Submission Portal. It describes experiment metadata.|
| Submitted File Metadata|submitted_file.tsv|Upload this file to the Submission Portal. It describes file metadata.|
| Manifest|manifest.yml|Download this file from the Submission Portal. This manifest is used by the GDC DTT for the actual file upload.|
| FAST File|ExperimentFile.fastq|Upload this file to the GDC DTT along with the manifest.|



## Step1. Prepare Files

The [GDC Dictionary](../../Dictionary/viewer.md) describes the types of entities that can be uploaded to the GDC.

The user can go to the GDC Dictionary to download the template files to be used for the upload. The templates can be populated with data by the user and should result in a valid file (if validation rules detailed in the dictionary are met).

A template file describes an entity with the following information:

* __Type__: identification of the entity.
* __IDs__: Project ID and Submitter ID of the entity.
* __Links__: Submitter ID of the links to other entities.
* __Properties__: properties of the entity.



Example of a __demographic file__ that can be uploaded in TSV format:

```tsv
type	project_id	submitter_id	cases.submitter_id	ethnicity	gender	race	year_of_birth	year_of_deathdemographic	TCGA-DEV3	TCGA-DEV-3-CASE-000-D1	TCGA-DEV-3-CASE-000	hispanic or latino	male	white	1950	0demographic	TCGA-DEV3	TCGA-DEV-3-CASE-001-D1	TCGA-DEV-3-CASE-001	not reported	female	white	1956	0
```


Once the user has prepared their files (in TSV or JSON format), they can move on to the next step, uploading their data through the Upload Data Wizard.

## Step 2. Upload Data Wizard


The GDC Data Submission Portal is equipped with a wizard window to guide you through the upload and validation of data. There are three stages:

* __Upload Files__: Upload a file into the user's browser, at this point nothing is submitted to the project workspace.
* __Validate Files__: Send a file to the GDC backend to validate its content (see below).
* __Confirm Submission__: Submit a validated file to the project workspace and produce a report.

The _'File Validation'_ stage acts as a safeguard against submitting incorrect files to the GDC Data Submission Portal. During the validation stage, the GDC API will validate the content of submitted files against the project's dictionary to detect potential errors. Invalid files will be flagged and the upload to the GDC will be denied until corrections are made by the user. A validation error report provided by the system can be used to isolate and correct errors for resubmission.

### Upload Files

From the project dashboard, clicking on _'UPLOAD'_ will open the submission wizard.

[![GDC Submission Wizard Upload Files](images/GDC_Submission_Wizard_Upload.png)](images/GDC_Submission_Wizard_Upload.png "Click to see the full image.")

Files can be added either by clicking on _'CHOOSE FILE(S)'_ or by using drag and drop.

### Validate Files

As soon as the first file is added, the wizard will move to the _'VALIDATE'_ section and the user can continue to add files.

[![GDC Submission Wizard Validate Files](images/GDC_Submission_Wizard_Validate.png)](images/GDC_Submission_Wizard_Validate.png "Click to see the full image.")

Once all files have been added, clicking on _'VALIDATE'_ will check if the files are valid for submission.

[![Confirm Submission](images/GDC_Submission_Wizard_Confirm.png)](images/GDC_Submission_Wizard_Confirm.png "Click to see the full image.")

If the upload contains invalid files, the user will not be able to submit the data and those files will need to be either corrected and re-uploaded or removed from the submission.

[![Invalid Files in a Submission](images/GDC_Submission_wizard_Invalid_Files.png)](images/GDC_Submission_wizard_Invalid_Files.png "Click to see the full image.")

Files can be removed from the submission by clicking on the _'garbage'_ icon related to the file.

Click on _'View Report'_ to see the error report which provides details about why files failed validation. This report is detailed in the [Data Submission Reports](Reports.md) section of the documentation.



### Confirm Upload

Once all files are valid for upload, clicking on _'Confirm Upload'_ will upload your files to the GDC project workspace.

**Note:** Uploaded files are not released immediately. The user must submit the data to the GDC and release the project to make data available on the GDC Data Portal.

[![Successful Submission](images/GDC_Submission_wizard_successful_submission.png)](images/GDC_Submission_wizard_successful_submission.png "Click to see the full image.")

You can then click on _'CLICK HERE TO VIEW THE TRANSACTION'_ to be redirected to the transaction list page.

## Step 3. GDC Data Transfer Tool


**Step 3 is applicable to Experiment Data only.**

The GDC Data Transfer Tool is used to upload the actual file.

Once the user has uploaded the metadata through the Upload Data Wizard (e.g. read_group and submitted file metadata), they will be able to download the manifest from the transaction report.

**Note:** You can also download the manifest from the Browse menu.

[![Transaction Manifest](images/GDC_Submission_Transactions_Get_Manifest.png)](images/GDC_Submission_Transactions_Get_Manifest.png "Click to see the full image.")

Use this manifest to upload your actual files to the GDC Data Transfer Tool. Please refer to the [GDC Data Transfer Tool's User Guide] (../Data_Transfer_Tool/Users_Guide/Getting_Started.md) for more information.

## Download previously submitted files

The transaction page, accessible through the Browse menu, lists all previous transactions in the project. The user can download submitted files in the details section of the screen by selecting a particular transaction.

[![Transaction Original Files](images/GDC_Submission_Transactions_Original_Files.png)](images/GDC_Submission_Transactions_Original_Files.png "Click to see the full image.")
