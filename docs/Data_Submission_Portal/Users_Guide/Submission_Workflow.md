# Submission Workflow

## Overview

The workflow diagram available below presents the main features available within GDC Submission Portal in order to submit data and release it to the research community on the GDC Data Portal.

[![GDC Data Submission Portal Workflow](images/GDC_Submission_Portal_Workflow.png)](images/GDC_Submission_Portal_Workflow.png "Click to see the full image.")


## Project Life Cyle

ADD DIAGRAM HERE

Summary:

* I can upload and validate Data if my project __is not__ REVIEW.
* I can review my project if it is OPEN. It will prevent users from uploading new data during the review period.
* I can submit data to the GDC if my project is REVIEW.
* I can release my project at any time: it will release the harmonized data to the GDC Data Portal.

## Upload, Submit and Release

### Upload and Validate Data
The Submitter will upload data to the project workspace and validate data with the GDC dictionary. At this point, data is not submitted to GDC yet.
[![GDC Data Submission Portal Workflow Upload](images/GDC_Submission_Portal_Workflow_Upload.png)](images/GDC_Submission_Portal_Workflow_Upload.png "Click to see the full image.")


### Review and Submit
When data in the project workspace is ready for processing, the Submitter or Project Owner has to submit the data to the GDC. It will trigger the [GDC Data Harmonization Pipeline](https://gdc.nci.nih.gov/submit-data/gdc-data-processing-software-and-algorithms/2-data-harmonization).

Two main actions should be performed:

* REVIEW the project: it will prevent other users from uploading new data to the project. The user should verify that data is ready for processing.
* SUBMIT data to the GDC: After reviewing the project data, the user can submit it to the GDC. It will trigger the harmonization process.

However, during the REVIEW process, if the user thinks the data is not ready for processing then he can RE-OPEN the project. He would then be able to upload more data the project workspace.


[![GDC Data Submission Portal Workflow Submit](images/GDC_Submission_Portal_Workflow_Submit.png)](images/GDC_Submission_Portal_Workflow_Submit.png "Click to see the full image.")

### Release
When GDC harmonized data is ready and project data is complete, the Project Owner will release the project. It will release harmonized data to the GDC Data Portal.


[![GDC Data Submission Portal Workflow Release](images/GDC_Submission_Portal_Workflow_Release.png)](images/GDC_Submission_Portal_Workflow_Release.png "Click to see the full image.")

