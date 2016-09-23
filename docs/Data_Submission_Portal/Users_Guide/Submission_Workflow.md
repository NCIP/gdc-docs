# Submission Workflow

## Overview

The workflow diagram below represents the data submission process that is implemented by the GDC Data Submission Portal. The submitter logs into the GDC Data Submission Portal, uploads data into the project workspace, and validates the data. When the data is ready for processing, the submitter reviews the data and submits it to the GDC. The GDC processes the data through the [GDC Harmonization Process](https://gdc.nci.nih.gov/submit-data/gdc-data-harmonization). Processed data is released to the [GDC Data Portal](https://gdc-portal.nci.nih.gov/) and other [GDC Data Access Tools](https://gdc.nci.nih.gov/access-data/data-access-processes-and-tools) according to [GDC data sharing policies](https://gdc.cancer.gov/submit-data/data-submission-policies).

[![GDC Data Submission Portal Workflow](images/gdc-submission-portal-submission-workflow.png)](images/gdc-submission-portal-submission-workflow.png "Click to see the full image.")

### Upload and Validate Data

The submitter uploads clinical and biospecimen data to the project workspace using GDC templates that are available in the [GDC Data Dictionary](https://gdc-docs.nci.nih.gov/Data_Dictionary/). The GDC will validate the uploaded data against the GDC Data Dictionary.

To upload data files such as BAM or FASTQ, the submitter registers file metadata with the GDC in a manner similar to uploading clinical and biospecimen data. After registering the files, the submitter downloads a manifest from the GDC Data Submission Portal and uses it with the [GDC Data Transfer Tool](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool) to upload the files. Alternatively, registered files can be uploaded directly to the [GDC API](https://gdc.cancer.gov/developers/gdc-application-programming-interface-api) without the use of the Data Transfer Tool.

[![GDC Data Submission Portal Workflow Upload](images/gdc-submission-portal-data-upload-workflow.png)](images/gdc-submission-portal-data-upload-workflow.png "Click to see the full image.")

### Review and Submit Data

When all necessary data and files have been uploaded, the submitter reviews the dataset and submits it to the GDC for processing through the [GDC Data Harmonization Pipeline](https://gdc.nci.nih.gov/submit-data/gdc-data-harmonization).

### Release Data

The GDC will release data according to [GDC data sharing policies](https://gdc.cancer.gov/submit-data/data-submission-policies). Data may be released after six months from the date of upload, or the submitter may request earlier release using the "Release Project" function.

Harmonized data will be released to GDC users through the [GDC Data Portal](https://gdc-portal.nci.nih.gov/) and other [GDC Data Access Tools](https://gdc.nci.nih.gov/access-data/data-access-processes-and-tools).

## Project Lifecycle

The lifecycle of a project in the GDC describes the workflow throughout the data submission process. The project lifecycle starts with the upload and validation of data into the project and ends with the release of the harmonized data to the GDC Data Portal and other GDC data access tools. Throughout the lifecycle, the project transitions through various states in which the project is open for uploading data, in review, and processing. This lifecycle is continuous as new project data becomes available.

To summarize the project lifecycle and transitions to the various states of the project, the following operations can be performed by the submitter:

*Note:* "Submit to the GDC" and "Release" actions can be performed only if the user has release privileges.

### Project State
The diagram below demonstrates the transition of a project through the various states. Initially the project is OPEN for data upload and validation. When the data is uploaded and ready for review, the submitter changes the project state to REVIEW. During the REVIEW state, the project is locked so that additional data cannot be uploaded. If data changes are needed during the review period, the project can be re-opened and the state changes back to OPEN. When review has been completed and the submitter submits the data for GDC processing, the project state changes to SUBMITTED. When the data has been processed, the project state changes back to OPEN for new data to be submitted to the project.

[![GDC Data Submission Portal Project State](images/gdc-submission-portal-project-states.png)](images/gdc-submission-portal-project-states.png "Click to see the full image.")


## File Status Lifecycle

This section describes states pertaining to experimental data files throughout the data submission process. An experimental data file could be genomic sequence file (such as a BAM or FASTQ) or a pathology slide image file. The file lifecycle starts when a submitter uploads metadata for a file to the GDC Data Submission Portal. Upload of file metadata registers a description of the file in the GDC. The submitter can then use the [GDC Data Transfer Tool](https://gdc.nci.nih.gov/access-data/gdc-data-transfer-tool) to upload the actual file. Throughout the lifecycle, the file status transitions through various states from when the file is initially registered through file submission and processing.

[![GDC Data Submission Portal File Status](images/gdc-submission-portal-file-state-vs-state.png)](images/gdc-submission-portal-file-state-vs-state.png "Click to see the full image.")
