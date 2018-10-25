# Getting Started

## Overview

The National Cancer Institute (NCI) Genomic Data Commons (GDC) Data Submission Portal User's Guide is the companion documentation for the [GDC Data Submission Portal](https://gdc.cancer.gov/submit-data/gdc-data-submission-portal) and provides detailed information and instructions for its use.

The GDC Data Submission Portal is a platform that allows researchers to submit and release data to the GDC. The key features of the GDC Data Submission Portal are:

* __Upload and Validate Data__: Project data can be uploaded to the GDC project workspace. The GDC will validate the data against the [GDC Data Dictionary](https://gdc-docs.nci.nih.gov/Data_Dictionary/).
* __Review and Submit Data__: Prior to submission, data can be reviewed to check for accuracy. Once the review is complete, the data can be submitted to the GDC for processing through [Data Harmonization](https://gdc.cancer.gov/submit-data/gdc-data-harmonization).
* __Release Data__: After harmonization, data can be released to the research community for access through [GDC Data Access Tools](https://gdc.cancer.gov/access-data/data-access-processes-and-tools).
* __Browse and Download Data__: Data that has been uploaded into the project workspace can be browsed to ensure that the data is ready for processing and downloaded for review or update. Data can then be re-uploaded before it is released for access through [GDC Data Access Tools](https://gdc.cancer.gov/access-data/data-access-processes-and-tools).
* __Status and Alerts__: Visual cues are implemented to easily identify incomplete submissions.
* __Transactions__: List of all actions performed in a project are displayed along with details about each transaction.

## Key Features

### Upload and Validate Data
To submit data to the GDC, the user will prepare the data and upload it to the project workspace.

The main categories of data that can be uploaded include:

* __Clinical Data__: Elements such as `gender`, `age`, `diagnosis`, etc. as defined in the GDC Data Dictionary.
* __Biospecimen Data__: Information about entities such as `samples`, `aliquots`, etc. as defined in the GDC Data Dictionary.
* __Submittable Data Files__: Sequencing data such as BAM and FASTQ files, slide images, and other experimental data collected by the study.

The [GDC Data Dictionary Viewer](../../Data_Dictionary/viewer.md) outlines the minimum field requirements for each of the three categories listed above.

### Review and Submit Data

Once data is uploaded to the project workspace, it can be reviewed to ensure that the data is ready for processing through the [GDC Harmonization Process](https://gdc.cancer.gov/submit-data/gdc-data-harmonization). The review will lock the project to ensure that additional data cannot be uploaded while in review. During this period the data can be browsed or downloaded in the Data Submission Portal.

If the project is ready for processing, data can be submitted to the GDC. If the project is not ready for processing, the project can be re-opened. This will allow for additional data to be uploaded to the project workspace.

### Release Data

The GDC will release data according to [GDC data sharing policies](https://gdc.cancer.gov/submit-data/data-submission-policies). Data must be released within six months after GDC data processing has been completed, or the submitter may request earlier release using the "Request Release" function.

Upon release, harmonized data will be available to GDC users through the [GDC Data Portal](https://portal.gdc.cancer.gov/) and other [GDC Data Access Tools](https://gdc.cancer.gov/access-data/data-access-processes-and-tools).

After harmonized data is released, it can only be redacted by GDC administrators under certain conditions. To request redaction of released data, please contact [GDC User Services](https://gdc.cancer.gov/support#gdc-help-desk).

### Browse and Download Data

Authorized submitters can browse and retrieve data submitted to their project using the Data Submission Portal.  Retrieval of data submitted to the submission portal can be accomplished by using the API or the Data Transfer Tool.  UUIDs of submitted files can be retrieved from the submission portal or with a [GraphQL](https://docs.gdc.cancer.gov/API/Users_Guide/Submission/#querying-submitted-data-using-graphql) query.   Please see the [API](https://docs.gdc.cancer.gov/API/Users_Guide/Downloading_Files/) documentation for more information about downloads.


### Status and Alerts

The GDC Data Submission Portal Dashboard and navigation panel displays a summary of submitted data and associated data elements, such as the number of cases with Clinical data or Biospecimen data.

### Transactions

Submitters can access a list of all actions performed in a project by clicking on the Transactions tab on the dashboard. This will display a list of all past transactions for the selected project. Users can access details about each transaction. The most recent transactions are also displayed on the dashboard.

## Release Notes

The [Release Notes](../../Data_Submission_Portal/Release_Notes/Data_Submission_Portal_Release_Notes.md) section of this User's Guide contains details about new features, bug fixes, and known issues.
