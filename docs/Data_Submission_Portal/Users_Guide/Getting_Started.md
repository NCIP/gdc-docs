# Getting Started

## Overview

The National Cancer Institute (NCI) Genomic Data Commons (GDC) Data Submission Portal User's Guide is the companion documentation for the [GDC Data Submission Portal](https://gdc.nci.nih.gov/submit-data/gdc-data-submission-portal) and provides detailed information and instructions for its use.

The GDC Data Submission Portal is a platform that allows researchers to submit and release data into the GDC. The key features of the GDC Data Submission Portal are:

* __Upload and Validate Data__: Upload project data to the GDC project workspace. The GDC will validate the data against the [GDC Data Dictionary](https://gdc-docs.nci.nih.gov/Data_Dictionary/).
* __Review and Submit Data__: Review project data which will lock the project to ensure that additional data cannot be uploaded while in review. Once the review is complete, the data can be submitted to the GDC for processing through the [GDC Harmonization Process](https://gdc.nci.nih.gov/submit-data/gdc-data-harmonization).
* __Release Data__: Release data to the research community for access through [GDC Data Access Tools](https://gdc.nci.nih.gov/access-data/data-access-processes-and-tools).
* __Download Data__: Download data that has been uploaded into the project workspace for review or update. Re-upload data before the data is released for access through [GDC Data Access Tools](https://gdc.nci.nih.gov/access-data/data-access-processes-and-tools).
* __Browse Data__: Browse data that has been uploaded to the project workspace to ensure that the project is ready for processing.
* __Status and Alerts__: Obtain visual cues to easily identify incomplete submissions.
* __Reports__: Access and download submission-related reports to gain better insight into the status of submitted data.

## Submission Workflow

The Submission Workflow is described in detail in the [Data Submission Workflow](Submission_Workflow.md) section. 

## Key Features

### Upload and Validate Data
In order to submit data to the GDC, the submitter will prepare the data and upload it to the project workspace.

The main categories of data that can be uploaded include: 

* Clinical Data: Elements such as gender, age, diagnosis, etc. as defined in the GDC Data Dictionary.
* Biospecimen Data: Entities such as samples, aliquots, etc. as defined in the GDC Data Dictionary.
* Experiment Data: At least one type of experiment data, such as a read group, as defined in the GDC Data Dictionary.

By using the [GDC Data Dictionary Viewer](../../Data_Dictionary/viewer.md), the user can identify the minimum field requirements for each of the three categories listed above.

In general, submitters have six months to upload data to the project workspace before it should be submitted to the GDC.

### Review and Submit Data

Once data is uploaded to the project workspace, the submitter or project owner should review it to ensure that the data is ready for processing through the ([GDC Harmonization Process](https://gdc.nci.nih.gov/submit-data/gdc-data-harmonization). The review will lock the project to ensure that additional data cannot be uploaded while in review. During this period, the submitter or project owner can browse or download the data in the Data Submission Portal. 

If the project is ready for processing, the submitter or project owner can submit the data to the GDC. If the project is not ready for processing, the submitter or project owner can re-open it. This will allow the submitter to upload more data to the project workspace.

In general, submitters are requested to submit their data to the GDC within six months from the first upload to the project workspace.

### Release Data

To ensure that the GDC always releases high-quality data, harmonized data is not automatically released. The Submitter or Project Owner will be able to review harmonized data and decide whether to release the data to the GDC Data Portal and other GDC data acess tools. The confirmation of the release will initiate data indexing and will release the data to the [GDC Data Portal](https://gdc-portal.nci.nih.gov/projects/t) and other GDC data access tools.

In general, submitters are requested to release files within six months from submission to the GDC. This guideline is in place in support of the GDC's aim to make data available to the community according to the [NCI Genomic Data Sharing Policy](http://www.cancer.gov/grants-training/grants-management/nci-policies/genomic-data).

Released cases and/or files can be redacted from the GDC. Redaction is performed by GDC administrators, at the case level through synchronization with dbGaP, and at file level through a submitter's request usually after a data quality issue is identified. The GDC Data Submission Portal itself currently does not support redaction through the web user interface.

### Browse and Download Data

In addition to submitting and releasing files, the GDC Data Submission Portal can also display submitted files and entities through multiple pages and tables.

Along with listing entities, the GDC Data Submission Portal also offers users the ability to download files submitted initially or the latest version of a specific entity (e.g. the clinical data for a case).

### Status and Alerts

The GDC Data Submission Portal dashboard and navigation panel displays the number of missing elements to be attached to an entity, such as cases missing clinical data or biospecimen data.

### Transactions

The user can access a list of all actions performed in a project by clicking on transactions in the Browse menu (e.g., upload data, review project, etc.). This will display a list of all past transactions for the selected project. By clicking on a transaction, users can access details of this transaction as well as download the original uploaded files.

Transactions are also displayed on the dashboard to highlight the most recent transactions.

### Reports

The GDC Data Submission Portal provides a broad range of submission oriented reports. Reports are detailed in the [Reports](Reports.md) section.

## Access

The GDC Data Submission Portal is accessible using a web browser such as Chrome, Internet Explorer, or Firefox at the following URL: [https://gdc-portal.nci.nih.gov/submission/](https://gdc-portal.nci.nih.gov/submission/).
Upon navigating to the site, the GDC Data Submission Portal login page or Dashboard is displayed.

## Release Notes

The GDC Data Submission Portal is regularly being updated with new features. The [Release Notes](../Release_Notes/Data_Submission_Portal_Release_Notes.md) section contains details about new features, bug fixes and known issues.
