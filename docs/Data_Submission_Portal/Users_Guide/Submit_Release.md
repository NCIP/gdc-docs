# Submit and Release

## Overview

The GDC Data Submission process is detailed on the [GDC Website]( https://gdc.nci.nih.gov/submit-data/data-submission-processes-and-tools).

This chapter will focus on the following features:

* The review and submission of data to the GDC.
* The Release of data to the GDC Data Portal

## Review and Submit

Once data is uploaded to the project workspace (see previous section: [Upload Data](Upload_Data.md), the submitter or project owner should review it to ensure that data is ready for processing by the GDC ([Harmonization Process](https://gdc.nci.nih.gov/submit-data/gdc-data-harmonization)).

The review will lock the project and ensure no user can upload more data. During that period, the submitter or project owner can browse the data in the Submission Portal or download it. 

If the project is ready for processing, the submitter or project owner will submit data to the GDC.

If the project is not ready for processing, the submitter or project owner can re-open it. Then the submitter will be able to upload more data to the project workspace.

In general, submitters are requested to submit their data to the GDC within six months from the first upload to the project workspace.

[![GDC Data Submission Portal Workflow Submit](images/GDC_Submission_Portal_Workflow_Submit.png)](images/GDC_Submission_Portal_Workflow_Submit.png "Click to see the full image.")

### Review

To review and submit data to the GDC, the user should have release privileges. The user will be able to view the section below on the dashboard.

The "Review" button is available only if the project is in "OPEN" state.

[![GDC Submission Review Tab](images/GDC_Submission_Submit_Release_Review_tab.png)](images/GDC_Submission_Submit_Release_Review_tab.png "Click to see the full image.")

Reviewing the project will prevent other users from uploading data to the project. Once the review is complete, the user can submit data to the GDC.

Once the user clicks on "Review", the project state will change to "REVIEW":

[![GDC Submission Review State](images/GDC_Submission_Submit_Release_Project_State_Review.png)](images/GDC_Submission_Submit_Release_Project_State_Review.png "Click to see the full image.")



### Submit to the GDC

The "Submit" button is available only if the project is in "REVIEW" state. At this point, the user can decide whether to re-open the project to upload more data or to submit data to the GDC.

Once the user submits data to the GDC, they __cannot upload additional data until the harmonization process is complete__.

[![GDC Submission Submit Tab](images/GDC_Submission_Submit_Release_Submit_tab.png)](images/GDC_Submission_Submit_Release_Submit_tab.png "Click to see the full image.")

When the user clicks on the action "Submit", the following Submit popup is displayed:

[![GDC Submission Submit Popup](images/GDC_Submission_Submit_Release_Submit_Popup.png)](images/GDC_Submission_Submit_Release_Submit_Popup.png "Click to see the full image.")


After the user clicks on "Submit Validated Data to the GDC", the project state becomes "Submitted":
[![GDC Submission Project State](images/GDC_Submission_Submit_Release_Project_State.png)](images/GDC_Submission_Submit_Release_Project_State.png "Click to see the full image.")


## Release
The project owner should release data when harmonization process is complete and the harmonized data is ready to release.
Releasing the data will make the data available to the GDC Data Portal.

As soon as the project owner releases the project, any harmonized (updated or added to the project) will be released to the [GDC Data Portal](https://gdc-portal.nci.nih.gov/projects/t).

[![GDC Data Submission Portal Workflow Release](images/GDC_Submission_Portal_Workflow_Release.png)](images/GDC_Submission_Portal_Workflow_Release.png "Click to see the full image.")

To release data to the GDC Data Portal, the user should have release privileges. The user will be able to view the following section on the dashboard.

[![GDC Submission Release Tab](images/GDC_Submission_Submit_Release_Release_tab.png)](images/GDC_Submission_Submit_Release_Release_tab.png "Click to see the full image.")

When the user clicks on the action "Release", the following Release popup is displayed:

[![GDC Submission Release Popup](images/GDC_Submission_Submit_Release_Release_Popup.png)](images/GDC_Submission_Submit_Release_Release_Popup.png "Click to see the full image.")

After the user clicks on "Release Submitted and Processed Data", the project release state becomes "Released":

[![GDC Submission Project State](images/GDC_Submission_Submit_Release_Project_State.png)](images/GDC_Submission_Submit_Release_Project_State.png "Click to see the full image.")


__Note__: Released cases and/or files can be redacted from the GDC. Redaction is performed by GDC administrators, at case level through synchronization with dbGaP, and at file level at the submitter's request usually after a data quality issue is identified. The GDC Data Submission Portal itself currently does not support redaction through the web user interface.
