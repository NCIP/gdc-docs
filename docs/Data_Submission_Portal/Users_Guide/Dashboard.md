# Dashboard

## Overview

The GDC Data Submission Portal dashboard provides details about a specific project.

[![GDC Submission Dashboard Page](images/GDC_Submission_Dashboard.png)](images/GDC_Submission_Dashboard.png "Click to see the full image.")

The dashboard contains various visual elements to guide the user through all stages of submission, from viewing the [Data Dictionary](https://gdc-docs.nci.nih.gov/Data_Dictionary/) in support of data upload and validation to releasing a project.

To better understand the data displayed on the dashboard and the available actions, please refer to the [Submission Workflow](Submission_Workflow.md).

## Project Status

The top section of the dashboard is broken down into four primary charts:

* __Cases with Clinical__: Details the number of cases for which clinical data have been uploaded.
* __Cases with Biospecimen__: Details the number of cases for which biospecimen data have been uploaded.
* __Cases with Experiment Data__: Details the number of cases for which experiment data have been uploaded.
* __Files Uploaded__: Details the number of files uploaded through the GDC Data Transfer Tool. For more information on this chart, please refer to [File Status Life Cycle](Submission_Workflow.md#file-status-life-cycle). 

These charts are constantly updated to reflect the current state of the selected project.

Clicking on "MORE" opens a table view with additional details.

[![GDC Submission Dashboard Details Widget](images/GDC_Submission_Dashboard_Details.png)](images/GDC_Submission_Dashboard_Details.png "Click to see the full image.")

## Action Tabs

There are three action tabs available in the middle section of the dashboard.

* [Upload & Validate](Upload_Data.md): Allows a submitter to upload project data to the GDC project workspace. The GDC will validate the uploaded data against the [GDC Data Dictionary](https://gdc-docs.nci.nih.gov/Data_Dictionary/).
* [Review & Submit](Submit_Release.md#review-and-submit): Allows a submitter to review project data which will lock the project to ensure that additional data cannot be uploaded while in review. Once the review is complete, the data can be submitted to the GDC for processing through the [GDC Harmonization Process](https://gdc.nci.nih.gov/submit-data/gdc-data-harmonization).
* [Release](Submit_Release.md#release): Allows a submitter to release data to the research community for access through [GDC Data Access Tools](https://gdc.nci.nih.gov/access-data/data-access-processes-and-tools).

These actions and associated features are further detailed in the respective sections of the documentation.

## Latest Transactions

The latest transactions section in the homepage lists the most recent [transactions](Transactions.md) associated with the projects that the user has access to.

## Reports

The reports section in the dashboard provides access to project reports on data submission, based on user authorization. More details about the reports are available in the Reports section.
