# Getting Started

## Overview

The National Cancer Institute (NCI) Genomic Data Commons (GDC) Data Submission Portal User's Guide is the companion documentation for the [GDC Data Submission Portal](https://gdc.nci.nih.gov/submit-data/gdc-data-submission-portal) and provides detailed information and instructions for its use.

The GDC Data Submission Portal is a platform allowing researchers to submit and release their data into GDC. The key features of the GDC Data Submission Portal are:

* __Upload and Validate Data__: Upload project data, metadata and associated files to the GDC project workspace. GDC will validate the data with the GDC data dictionary.
* __Review and Submit__: Review your project data will lock the project and ensure no user can upload more data. Once project data is reviewed, Submit it to GDC for data processing (Harmonization Pipeline).
* __Release__: Release your harmonized data to the research community (on the GDC Data Portal).
* __Download Data__: Download submitted data from your project workspace to review them or update them for re-upload before their release to the GDC Data Portal.
* __Browse Data__: Browse submitted data in the GDC Data Submission Portal to make sure your project is ready for processing.
* __Status and Alerts__: Visual mechanisms to easily identify incomplete submissions.
* __Reports__: Access and download submission-related reports to get better insight into the status of submitted data.

## Submission Workflow

The Submission Workflow is described in this section: [Data Submission Workflow](Submission_Workflow.md). 


## Key Features

### Upload and Validate
TO BE COMPLETED

### Review and Submit
TO BE COMPLETED

### Release

To ensure GDC always releases high-quality data, submitted files are not automatically released.

Once cases meet a minimum set of requirements, a submitter can sign-off on those cases via the submission portal. The sign-off will initiate data indexing and harmonization in preparation for a release of the data to the [GDC Data Portal] (https://gdc-dev.nci.nih.gov/access-data/gdc-data-portal).

Minimum requirements for releasing a case vary from one project to another, but must include data from the following three categories:

* Clinical Data: elements such as gender, age, diagnosis, etc. as defined in the GDC dictionary.
* Biospecimen Data: entities such as sample, aliquot etc as defined in the GDC dictionary.
* Molecular Data: at least one type of molecular data, such as, lane level DNA Sequences, as defined in the GDC dictionary.

By using the dictionary viewer, the user can identify minimum fields requirements for each of the three categories listed above.

In general, submitters are requested to release their files within six months from first submission. This guideline is in place in support of GDC's aim to make data available to the community according to the Genomic Data Sharing policies outlined by NCI.

Released cases and/or files can be redacted from GDC. Redaction is performed by GDC administrators, at case level through synchronization with dbGaP, and at file level through submitter's request usually after a data quality issue is identified. The GDC Data Submission Portal itself currently does not support redaction by researchers through the web user interface.

### Data Access and Download

In addition to submitting and releasing files, the GDC Data Submission Portal can also display submitted files and entities through multiple pages and tables.

Along with listing entities, the GDC Data Submission Portal also offers users the ability to download files submitted initially or the latest version of a specific entity (such as case).

### Status and Alerts

Via its Dashboard and navigation panel, the GDC Data Submission Portal displays the number of missing elements to be attached to an entity, such as cases missing clinical data or aliquots missing data bundles.

Displaying all of this on the same page provides a snapshot of the overall submission status of one or more projects.

### Dictionary Viewer

Although all projects are attached to the same data model, custom and project-specific dictionaries are added to extend the data model to the specific vocabulary and data elements of a given project.

To provide users with a better understanding of dictionaries attached to projects, an online dictionary viewer is available within each project. This allows users to adapt their submission files to the project they are submitting to.


### Transactions

The user can access a list of all previous submissions by clicking on transactions in the left panel. This will display a list of all past transactions for the selected project. By clicking on a transaction users can access details of this transaction as well as download the submitted files.

Transactions are also displayed on the dashboard via a widget displaying the most recent transactions.

### Reports

The GDC Data Submission Portal provides a broad range of submission oriented reports. Those reports are detailed in the [Data Submission Reports](Reports.md) section of the documentation.

## Access

The GDC Data Submission Portal is accessible using a web browser such as Chrome, Internet Explorer or Firefox at the following URL: [https://gdc-portal.nci.nih.gov/submission/](https://gdc-portal.nci.nih.gov/submission/).
Upon loading the site, the GDC Data Submission Portal login page or Dashboard is displayed.

## Release Notes

The GDC Data Portal is regularly being updated with new features. The [Release Notes](../Release_Notes/Data_Submission_Portal_Release_Notes.md) sections of the documentation contains details about new features, bug fixes and known issues.
