# Getting Started

## Overview

The National Cancer Institute (NCI) Genomic Data Commons (GDC) Data Submission Portal User's Guide is the companion documentation for the [GDC Data Submission Portal](https://gdc.nci.nih.gov/submit-data/gdc-data-submission-portal) and provides detailed information and instructions for its use.

The GDC Data Submission Portal is a platform that allows researchers to submit and release their data into the GDC. The key features of the GDC Data Submission Portal are:

* __Upload and Validate Data__: Upload project data to the GDC project workspace. The GDC will validate the data with the GDC data dictionary.
* __Review and Submit__: Review your project data will lock the project and ensure no user can upload more data. Once project data is reviewed, Submit it to the GDC for data processing (Harmonization Process).
* __Release__: Release your harmonized data to the research community (on the GDC Data Portal).
* __Download Data__: Download uploaded data from your project workspace to review them or update them for re-upload before their release to the GDC Data Portal.
* __Browse Data__: Browse uploaded data in the GDC Data Submission Portal to make sure your project is ready for processing.
* __Status and Alerts__: Visual mechanisms to easily identify incomplete submissions.
* __Reports__: Access and download submission-related reports to gain better insight into the status of submitted data.

## Submission Workflow

The Submission Workflow is described in this section: [Data Submission Workflow](Submission_Workflow.md). 


## Key Features

### Upload and Validate
In order to submit data to the GDC, the submitter will prepare the data and upload it to the project workspace.

These are the main categories of data that can be uploaded: 

* Clinical Data: elements such as gender, age, diagnosis, etc. as defined in the GDC Dictionary.
* Biospecimen Data: entities such as sample, aliquot, etc. as defined in the GDC Dictionary.
* Experiment Data: at least one type of Experiment data, such as a Read Group, as defined in the GDC Dictionary.

By using the [GDC Dictionary Viewer](../../Dictionary/viewer.md), the user can identify minimum fields requirements for each of the three categories listed above.


In general, submitters have six month to upload data to the project workspace before they should submit it to the GDC (see next section).

### Review and Submit
Once data is uploaded to the project workspace, the submitter or project owner should review it to ensure that data is ready for processing by the GDC ([Harmonization process](https://gdc.nci.nih.gov/submit-data/gdc-data-processing-software-and-algorithms/2-data-harmonization)).

The review will lock the project and ensure no user can upload more data. During that period, the submitter or project owner can browse the data in the Submission Portal or download it. 

If the project is ready for processing, the submitter or project owner will submit data to the GDC.

If the project is not ready for processing, the submitter or project owner can re-open it. Then the submitter will be able to upload more data to the project workspace.

In general, submitters are requested to submit their data to the GDC within six months from the first upload to the project workspace.


### Release

To ensure the GDC always releases high-quality data, harmonized data is not automatically released.

The Submitter or Project Owner will be able to review harmonized data then can decide to release the data to the GDC Data Portal. The confirmation of the release will initiate data indexing and will release data to the [GDC Data Portal](https://gdc-portal.nci.nih.gov/projects/t).

In general, submitters are requested to release their files within six months from the submission to the GDC. 

This guideline is in place in support of the GDC's aim to make data available to the community according to the Genomic Data Sharing policies outlined by NCI.

Released cases and/or files can be redacted from the GDC. Redaction is performed by GDC administrators, at case level through synchronization with dbGaP, and at file level through submitter's request usually after a data quality issue is identified. The GDC Data Submission Portal itself currently does not support redaction through the web user interface.

### Data Access and Download

In addition to submitting and releasing files, the GDC Data Submission Portal can also display submitted files and entities through multiple pages and tables.

Along with listing entities, the GDC Data Submission Portal also offers users the ability to download files submitted initially or the latest version of a specific entity (such case's clinical data).

### Status and Alerts

Via its Dashboard and navigation panel, the GDC Data Submission Portal displays the number of missing elements to be attached to an entity, such as cases missing clinical data or biospecimen data.


### Transactions

The user can access a list of all actions performed in a project by clicking on transactions in the Browse menu (e.g. upload data, review project, etc.). This will display a list of all past transactions for the selected project. By clicking on a transaction users can access details of this transaction as well as download the original uploaded files.

Transactions are also displayed on the dashboard by a widget displaying the most recent transactions.

### Reports

The GDC Data Submission Portal provides a broad range of submission oriented reports. Those reports are detailed in the [Data Submission Reports](Reports.md) section of the documentation.

## Access

The GDC Data Submission Portal is accessible using a web browser such as Chrome, Internet Explorer or Firefox at the following URL: [https://gdc-portal.nci.nih.gov/submission/](https://gdc-portal.nci.nih.gov/submission/).
Upon loading the site, the GDC Data Submission Portal login page or Dashboard is displayed.

## Release Notes

The GDC Data Submission Portal is regularly being updated with new features. The [Release Notes](../Release_Notes/Data_Submission_Portal_Release_Notes.md) section of the documentation contains details about new features, bug fixes and known issues.
