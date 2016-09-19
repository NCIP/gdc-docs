# Dashboard

## Overview

The GDC Data Submission Portal dashboard provides details about a specific project.

[![GDC Submission Dashboard Page](images/GDC_Submission_Dashboard.png)](images/GDC_Submission_Dashboard.png "Click to see the full image.")

The dashboard contains various visual elements to guide the user through all stages of submission, from viewing the [Data Dictionary](https://gdc-docs.nci.nih.gov/Data_Dictionary/) in support of data upload and validation to releasing a project.

To better understand the data displayed on the dashboard and the available actions, please refer to the [Submission Workflow](Submission_Workflow.md).

## Project Overview
The Project Overview sections of the dashboard displays the project state (open/ review/ submitted/ processing) and the GDC Release, which is XXXXXXXXXXXXXXX.

The remaining part of the top section of the dashboard is broken down into four primary charts:

* __Cases with Clinical__: Details the number of cases for which clinical data has been uploaded.
* __Cases with Biospecimen__: Details the number of cases for which biospecimen data has been uploaded.
* __Cases with Submittable Data Files__: Details the number of cases for which experiment data has been uploaded.
* __Submittable Data Files__: Details the number of files uploaded through the GDC Data Transfer Tool. For more information on this chart, please refer to [File Status Lifecycle](Submission_Workflow.md#file-status-life-cycle). The "DOWNLOAD MANIFEST" button below this chart allows the users to download a manifest for registered files in this project that have not yet been uploaded.

These charts are constantly updated to reflect the current state of the selected project.




[![GDC Submission Dashboard Details Widget](images/GDC_Submission_Dashboard_Details.png)](images/GDC_Submission_Dashboard_Details.png "Click to see the full image.")

## Action Panels

There are two action panels available below the Project Overview.

* [Upload & Validate](Upload_Data.md): Allows a submitter to upload project data to the GDC project workspace. The GDC will validate the uploaded data against the [GDC Data Dictionary](https://gdc-docs.nci.nih.gov/Data_Dictionary/). This panel also displays a table that displays details about the five latest [transactions](Transactions.md).
* [Review & Submit](Submit_Release.md#review-and-submit): Allows a submitter to review project data which will lock the project to ensure that additional data cannot be uploaded while in review. Once the review is complete, the data can be submitted to the GDC for processing through the [GDC Harmonization Process](https://gdc.nci.nih.gov/submit-data/gdc-data-harmonization).

These actions and associated features are further detailed in the respective sections of the documentation.
