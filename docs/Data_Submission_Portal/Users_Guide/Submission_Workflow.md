# Submission Workflow

## Overview

The workflow diagram below represents the data submission process implemented by the GDC Data Submission Portal. The submitter logs into the GDC Data Submission Portal, uploads data into the project workspace, and validates the data. When the data is ready for processing, the submitter reviews the data and requests submission for the data to GDC. The GDC harmonizes the data through the [GDC Harmonization Process](https://gdc.cancer.gov/submit-data/gdc-data-harmonization). The submitter can then release the harmonized data to the [GDC Data Portal](https://portal.gdc.cancer.gov/) and other [GDC Data Access Tools](https://gdc.cancer.gov/access-data/data-access-processes-and-tools) according to [GDC Data Sharing Policies](https://gdc.cancer.gov/submit-data/data-submission-policies).

[![GDC Data Submission Portal Workflow](images/gdc-submission-portal-submission-workflow.png)](images/gdc-submission-portal-submission-workflow.png "Click to see the full image.")

### Upload and Validate Data

The submitter uploads Clinical and Biospecimen data to the project workspace using GDC templates that are available in the [GDC Data Dictionary](https://docs.gdc.cancer.gov/Data_Dictionary/). The GDC will validate the uploaded data against the GDC Data Dictionary.

To upload submittable data files, such as sequence data in BAM or FASTQ format, the submitter must register file metadata with the GDC using a method similar to uploading Clinical and Biospecimen data. When files are registered, the submitter downloads a manifest from the GDC Data Submission Portal and uses it with the [GDC Data Transfer Tool](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool) to upload the files.

[![GDC Data Submission Portal Workflow Upload](images/gdc-submission-portal-data-upload-workflow.png)](images/gdc-submission-portal-data-upload-workflow.png "Click to see the full image.")

### Review and Submit Data

Once data is uploaded to the project workspace, it can be reviewed to ensure that the data is ready for processing through the [GDC Data Harmonization Pipeline](https://gdc.cancer.gov/submit-data/gdc-data-harmonization). Placing a project in REVIEW will close the project to ensure that additional data cannot be uploaded. During this period the data can be browsed or downloaded in the Data Submission Portal. If the project is ready for processing, data submission to the GDC can be requested. If the project is not ready for processing, the project can be re-opened. This will allow for additional data to be uploaded to the project workspace.



### Release Data

The GDC will release data according to [GDC data sharing policies](https://gdc.cancer.gov/submit-data/data-submission-policies). A project must be released within six months of harmonization or the submitter may request earlier release using the "Request Release Project" function.  

Harmonized data will be released to GDC users through the [GDC Data Portal](https://gdc-portal.nci.nih.gov/) and other [GDC Data Access Tools](https://gdc.cancer.gov/access-data/data-access-processes-and-tools). Once a project is released, all additional submitted data will automatically be released after harmonization.

## Project Lifecycle

The lifecycle of a project in the GDC describes the workflow throughout the data submission process. The project lifecycle starts with the upload and validation of data into the project and ends with the release of the harmonized data to the GDC Data Portal and other GDC data access tools. Throughout the lifecycle, the project transitions between being open for uploading data or in review where the project is closed and data cannot be uploaded. This lifecycle is continuous as new project data becomes available.

The diagram below demonstrates the lifecycle of a project. Initially the project is open for data upload and validation. Any changes to the data must be made while the project status is OPEN. When the data is uploaded and ready for review, the submitter changes the project status to REVIEW by clicking the REVIEW button on the dashboard. While in REVIEW, the project is closed so that additional data cannot be uploaded. If data changes are needed during the review period, the project can be re-opened and the and the status changes back to OPEN. When review has been completed and the submitter requests data submission to the GDC for processing. The project can be re-opened for new data to be submitted to the project.

[![GDC Data Submission Portal Project State](images/gdc-submission-portal-project-states.png)](images/gdc-submission-portal-project-states.png "Click to see the full image.")


## File Status Lifecycle

This section describes status pertaining to submittable data files throughout the data submission process. A submittable data file could contain data such as genomic sequences (such as a BAM or FASTQ) or pathology slide images. The file lifecycle starts when a submitter uploads metadata for a file to the GDC Data Submission Portal. Upload of file metadata, as a GDC Data Model entity, registers a description of the file in the GDC. The submitter can then use the [GDC Data Transfer Tool](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool) to upload the actual file. Throughout the lifecycle, the file status transitions through various states from when the file is initially registered through file submission and release. The diagram below details these status transitions.

[![GDC Data Submission Portal File Status](images/gdc-submission-portal-file-state-vs-state.png)](images/gdc-submission-portal-file-state-vs-state.png "Click to see the full image.")
