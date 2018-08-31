# Before Submitting Data to the GDC Portal

## Overview
The National Cancer Institute (NCI) Genomic Data Commons (GDC) Data Submission Portal User's Guide is the companion documentation for the [GDC Data Submission Portal](https://gdc.cancer.gov/submit-data/gdc-data-submission-portal) and provides detailed information and instructions for its use.

## Steps to Submit Data to the GDC
The following tasks are required to submit data to the [GDC Data Portal](https://portal.gdc.cancer.gov/).


* __Step 1:__ Complete the GDC Data [Submission Request Form](https://gdc.cancer.gov/data-submission-request-form). After submission, the requst will be reviewed by the GDC. 

   * __Step 1a:__ Durring this time create an eRA Commons account if you do not already have one. Follow the instructions found at the [NIH eRA Commons](https://era.nih.gov/registration_accounts.cfm).


* __Step 2:__ If the study is approved, contact a [Genomic Program Administrator (GPA)](https://osp.od.nih.gov/genomic-program-administrators/) to register the approved study in dbGaP.

   * __Step 2a__: Work with the GPA to register the study at the NCBI [database of Genotypes and Phenotypes (dbGaP)](https://www.ncbi.nlm.nih.gov/sra/docs/submitdbgap) as a trusted partner. When creating a project, give the project a two part identifiier, where the first portion is the __Program__ followed by a hyphen (__-__) and the second portion is the __Project__. 

   >*Example: The project __TCGA-BRCA__, the __Program__ is The Cancer Genome Atlas (__TCGA__), and their __Project__ is Breast Invasive Carcinoma (__BRCA__)*.

   * __Step 2b__: Provide the GPA with the approved final version of the [Institutional Certification](https://osp.od.nih.gov/scientific-sharing/institutional-certifications/), this will be passed onto dbGaP.

   * __Step 3b__: You may also provide the GPA with a list of Approved Data Submitters with eRA Common usernames.


* __Step 3:__ Accept the invitation to the dbGaP submission protal. Submit the study configuration file and subject IDs to the dbGaP study. When this is complete, the study will have a phs number.

* __Step 4:__ Before uploading data, familiarize yourself with the [data model](INSERT_URL) and [Data Dictionary](../../Data_Dictionary/viewer.md) to understand how different data relate to each other and where each data type is located.

* __Step 5:__ Go to the [GDC Data Portal](https://portal.gdc.cancer.gov/) and login by clicking the login button in the upper right hand corner. A pop up window will appear asking for eRA Commons Account credentials.

* __Step 6:__ Go to the [Data Submission Portal](https://portal.gdc.cancer.gov/submission/) to start uploading and validating data. For more information on this subject, go to [Data Submission Process](INSERT_URL).

## Key Features
The GDC Data Submission Portal is a platform that allows researchers to submit and release data to the GDC. The key features of the GDC Data Submission Portal are:

* __Upload and Validate Data__: Project data can be uploaded to the GDC project workspace. The GDC will validate the data against the [GDC Data Dictionary](../../Data_Dictionary/viewer.md).
* __Browse Data__: Data that has been uploaded to the project workspace can be browsed to ensure that the project is ready for processing.
* __Download Data__: Data that has been uploaded into the project workspace can be downloaded for review or update by using the [API](https://docs.gdc.cancer.gov/API/Users_Guide/Downloading_Files/) or the [Data Transfer Tool](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool). Data can then be re-uploaded before it is released for access through [GDC Data Access Tools](https://gdc.cancer.gov/access-data/data-access-processes-and-tools).
* __Review and Submit Data__: Prior to submission, data can be reviewed to check for accuracy. Once the review is complete, the data can be submitted to the GDC for processing through [Data Harmonization](https://gdc.cancer.gov/submit-data/gdc-data-harmonization).
* __Release Data__: After harmonization, data can be released to the research community for access through [GDC Data Portal](https://portal.gdc.cancer.gov/) and other [GDC Data Access Tools](https://gdc.cancer.gov/access-data/data-access-processes-and-tools).
* __Redaction__: After the release of harmonized data, if the need arises it can still be redacted by GDC administrators. To request redaction of released data, please contact [GDC User Services](https://gdc.cancer.gov/support#gdc-help-desk).
* __Status and Alerts__: Visual cues are implemented in the GDC Data Submission Portal Dashboard to easily identify incomplete submissions via panel displays summarizing submitted data and associated data elements.
* __Transactions__: A list of all actions performed in a project, with detailed information for each action.


## Submission Project Examples

Step-by-step instructions on GDC data submission and their relationship to the GDC Data Model are detailed in the [Data Submission Walkthrough](Data_Submission_Walkthrough.md).

## Release Notes

The [Release Notes](../../Data_Submission_Portal/Release_Notes/Data_Submission_Portal_Release_Notes.md) section of this User's Guide contains details about new features, bug fixes, and known issues.

