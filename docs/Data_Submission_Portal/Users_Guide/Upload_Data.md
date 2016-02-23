# Overview

GDC Data Submission process is detailed on the [GDC Website]( https://gdc.nci.nih.gov/submit-data/data-submission-processes-and-tools).
This chapter will focus on the upload and validation of the data.

# Types of files for upload
The GDC Data Submission Portal supports four types of files for submission:

* __Biospecimen__: Metadata describing a tissue specimen collected from a case and other material derived from samples for analysis.
* __Clinical__: A case’s clinical data.
* __Experiment Data__: GDC's unit of submission (see below).
* __Annotations__: Observations associated with any entity, which can be useful for interpreting the data.

More details about the submission process, data files and file formats can be found on the GDC website at [Data Submission Processes and Tools] (https://gdc.nci.nih.gov/submit-data/data-submission-processes-and-tools) and [Data Types and File Formats] (https://gdc.nci.nih.gov/submit-data/gdc-data-types-and-file-formats).

## Experiment Data

GDC has developed a submission unit called a “data bundle” (a set of files with associated metadata).

Data bundle types are defined by what files are expected and what introspection is done for validation or linking to other GDC entities. Each data bundle will be validated via a bundle type and project specific schema, including a JavaScript Object Notation (JSON) data dictionary, relationship check and molecular data quality check. 

The existing data bundles are described in the [GDC Dictionary] (../../Dictionary/viewer.md).


Following a data bundle upload, a __manifest__ can be obtained from the application and used to upload the actual files via the GDC Data Transfer Tool (DTT).

## Upload and Validation Process

Data upload and validation is one of the main feature of the GDC Data Submission Portal. Using a submission wizard, submitters are guided through a three-stage process:

* __File Upload__: Upload a file into the user's browser, at this point nothing is submitted to the project workspace.
* __File Validation__: Send the file to the GDC backend to validate its content (see below).
* __File Submission__: Submit validated file to the project workspace and produce a submission report.

The _'File Validation'_ stage acts as a safeguard against submitting incorrect files to the GDC Data Submission Portal. During the validation stage, the GDC API will validate content of submitted files against the project's dictionary to detect potential errors. Invalid files will be flagged and submission to GDC will be denied until corrections are made by the user. A validation error report provided by the system can be used to isolate and correct errors for resubmission.

The GDC Data Submission Portal supports the following file formats for submission:

* JSON
* TSV

During the upload and validation process, files are converted by the GDC API into entities and inserted into the database, maintaining a file-agnostic backend.

The GDC Data Submission Portal offers the ability to download files in different formats. To do so the system converts database entities back to the requested file format.


# Step1. Prepare files with the data dictionary
TO BE COMPLETED



# Step 2. Upload Data Wizard

TO BE UPDATED WITH NEW SCREENSHOTS


The GDC Data Submission Portal is equipped with a submission wizard window to guide you through the submission of the users files. When entering the submission process, this wizard window will provide an intuitive environment throughout the three stages of the submission:

* Upload Files.
* Validate Files.
* Submit Files.

## Upload Files

When in a project, clicking on _'SUBMIT'_ will open the submission wizard.

[![GDC Submission Wizard Upload Files](images/GDC_Submission_wizard_upload.png)](images/GDC_Submission_wizard_upload.png "Click to see the full image.")

Files can be added either by clicking on _'CHOOSE FILE(S)'_ or by using drag and drop.

## Validate Files

As soon as the first file is added, the wizard will move to the _'VALIDATE'_ section and the user can continue to add files.

[![GDC Submission Wizard Validate Files](images/GDC_Submission_wizard_Validate.png)](images/GDC_Submission_wizard_Validate.png "Click to see the full image.")

Once all files have been added, clicking on _'VALIDATE'_ will check if files are valid for submission.

[![Invalid Files in a Submission](images/GDC_Submission_wizard_invalid_files.png)](images/GDC_Submission_wizard_invalid_files.png "Click to see the full image.")

If the upload contains invalid files, the user will not be able to submit the data and those files will need to be either corrected and re-uploaded or removed from the submission.

Files can be removed from the submission by clicking on the _'garbage'_ icon related to the file.

An error report is also available to provide details about why files failed validation. This report is detailed in the [Data Submission Reports](https://gdc.nci.nih.gov/node/8449/) section of the documentation.

The dictionary viewer can be used to generate a template file to be used for submission. The template can be populated with data by the user and should result in a valid file (if validation rules detailed in the dictionary are met).

## Submit Files

Once all files are valid for submission, clicking on _'SUBMIT FILES'_ will submit the files to GDC.

**Note:** Submitted files are not released immediately. The user must select to release the files and GDC must validate and harmonize the submitted data for files to be released.

[![Successful Submission](images/GDC_Submission_wizard_successful_submission.png)](images/GDC_Submission_wizard_successful_submission.png "Click to see the full image.")

You can then click on _'CLICK HERE TO VIEW THE TRANSACTION'_ to be redirected to the transaction list page.

# Step 3. Upload Files through the DTT
TO BE COMPLETED


# Download previously submitted files

The transaction page, accessible through the left panel menu, list all previous transactions for a specified project. The user can download submitted files in the details section of the screen by selected a particular transaction.
