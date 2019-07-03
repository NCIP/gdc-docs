# Data Submission Portal Release Notes

| Version | Date |
|---|---|
| [v2.2.0](Data_Submission_Portal_Release_Notes.md#release-220) | February 20, 2019 |
| [v2.1.0](Data_Submission_Portal_Release_Notes.md#release-210) | November 7, 2018 |
| [v2.0.0](Data_Submission_Portal_Release_Notes.md#release-200) | August 23, 2018 |
| [v1.9.0](Data_Submission_Portal_Release_Notes.md#release-190) | May 21, 2018 |
| [v1.8.0](Data_Submission_Portal_Release_Notes.md#release-180) | February 15, 2018 |
| [v1.7.0](Data_Submission_Portal_Release_Notes.md#release-170) | November 16, 2017 |
| [v1.6.0](Data_Submission_Portal_Release_Notes.md#release-160) | August 22, 2017 |
| [v1.5.1](Data_Submission_Portal_Release_Notes.md#release-151) | March 16, 2017 |
| [v1.3.0](Data_Submission_Portal_Release_Notes.md#release-130) | October 31, 2016 |
| [v1.2.2](Data_Submission_Portal_Release_Notes.md#release-122) | September 23, 2016 |
| [v1.1.0](Data_Submission_Portal_Release_Notes.md#release-110) | May 20th, 2016 |
| [v0.3.24.1](Data_Submission_Portal_Release_Notes.md#release-03241) | February 26, 2016 |
| [v0.3.21](Data_Submission_Portal_Release_Notes.md#release-0321) | January 27, 2016 |
| [v0.2.18.3](Data_Submission_Portal_Release_Notes.md#release-02183) | November 30, 2015 |

## Release 2.2.0

* __GDC Product__: GDC Data Submission Portal
* __Release Date__: February 20, 2019

### New Features and Changes

* Renamed the "Request Submission" button to "Request Harmonization" to make the purpose of this action more clear. <!--SUBP-500-->

### Bugs Fixed Since Last Release

*  Fixed the right scroll bar in the records list on the Browse page so that it works in Firefox. <!--SUBP-493-->
*  Fixed a dead link to the Submission Portal User Guide on the Dashboard. <!--SUBP-494-->

### Known Issues and Workarounds

*  When creating entities in the Submission Portal, occasionally an extra transaction will appear with status error. This does not seem to impact that actual transaction, which is recorded as occurring successfully.
<!--API-219-->

## Release 2.1.0

* __GDC Product__: GDC Data Submission Portal
* __Release Date__: November 7, 2018

### New Features and Changes

* Updated the project columns to include a Release column in addition to the Batch Submit column. <!--SUBP-487-->

### Bugs Fixed Since Last Release

*  Fixed quick search so that projects with a dash in the name will no longer break the search. <!--SUBP-465-->
*  PO reports will now return the latest data for each project that has completed running. <!--TT-826-->

### Known Issues and Workarounds

*  When creating entities in the Submission Portal, occasionally an extra transaction will appear with status error. This does not seem to impact that actual transaction, which is recorded as occurring successfully.
<!--API-219-->

## Release 2.0.0

* __GDC Product__: GDC Data Submission Portal
* __Release Date__: August 23, 2018

### New Features and Changes

* The submission process has been updated to request submission and release. <!--SUBP-443-->
   * Users can "Request Submission" once data has been reviewed.  Previously the button that said "Submit" now says "Request Submission".
   * Users can "Request Release" once data has processed by the GDC.  Previously the button that said "Release" now says "Request Release".

### Bugs Fixed Since Last Release

*  Fixed Download All Manifest button being much larger than other buttons. <!--SUBP-461-->
*  Fixed bug where all tables said Showing x of xx projects instead of the correct entity. <!--SUBP-462-->

### Known Issues and Workarounds

*  When creating entities in the Submission Portal, occasionally an extra transaction will appear with status error. This does not seem to impact that actual transaction, which is recorded as occurring successfully.
<!--API-219-->


## Release 1.9.0

* __GDC Product__: GDC Data Submission Portal
* __Release Date__: May 21, 2018

### New Features and Changes

*  Added the ability to download all case metadata in the Browse Cases view <!--SUBP-443-->

### Bugs Fixed Since Last Release

*  Time outs when loading submission portal project list <!--SUBP-407-->
*  Missing PO reports for CPTAC-3 project <!--SUBP-407-->
*  Windows - Scroll bar interfering with browser scroll bar <!--SUBP-424-->
*  Donut "Cases with submittable data files" always shows 0 <!--SUBP-425-->
*  For Biospecimen entities, the "Download all" button does not take filtering into account <!--SUBP-423-->
*  For Clinical entities, the "Download all" button was downloading all clinical entities instead of selected clinical type <!--SUBP-442-->

### Known Issues and Workarounds

*  When creating entities in the Submission Portal, occasionally an extra transaction will appear with status error. This does not seem to impact that actual transaction, which is recorded as occurring successfully.
<!--API-219-->
Release details are maintained in the [GDC Data Submission Portal Change Log](https://github.com/NCI-GDC/submission-ui/blob/master/CHANGELOG.md).

## Release 1.8.0

* __GDC Product__: GDC Data Submission Portal
* __Release Date__: February 15, 2018

### New Features and Changes

*  Added the ability to support banners with hyperlinks <!--SUBP-411-->

### Bugs Fixed Since Last Release

*  Fixed 508 compliance issues <!--SUBP-408-->

### Known Issues and Workarounds

*  When creating entities in the Submission Portal, occasionally an extra transaction will appear with status error. This does not seem to impact that actual transaction, which is recorded as occurring successfully.
<!--API-219-->
Release details are maintained in the [GDC Data Submission Portal Change Log](https://github.com/NCI-GDC/submission-ui/blob/master/CHANGELOG.md).


## Release 1.7.0

* __GDC Product__: GDC Data Submission Portal
* __Release Date__: November 16, 2017

### New Features and Changes

*  None

### Bugs Fixed Since Last Release

*  Fixed bug where error would be produced even while project was successfully submitted <!--SUBP-394-->

### Known Issues and Workarounds

*  When creating entities in the Submission Portal, occasionally an extra transaction will appear with status error. This does not seem to impact that actual transaction, which is recorded as occurring successfully.
<!--API-219-->
Release details are maintained in the [GDC Data Submission Portal Change Log](https://github.com/NCI-GDC/submission-ui/blob/master/CHANGELOG.md).


## Release 1.6.0

* __GDC Product__: GDC Data Submission Portal
* __Release Date__: August 22, 2017

### New Features and Changes

*	 Added ability to see metadata for particular harmonized data files in the Submission Portal <!--SUBP-393-->

### Bugs Fixed Since Last Release

None

### Known Issues and Workarounds

*  When creating entities in the Submission Portal, occasionally an extra transaction will appear with status error. This does not seem to impact that actual transaction, which is recorded as occurring successfully.
<!--API-219-->
Release details are maintained in the [GDC Data Submission Portal Change Log](https://github.com/NCI-GDC/submission-ui/blob/master/CHANGELOG.md).


## Release 1.5.1

* __GDC Product__: GDC Data Submission Portal
* __Release Date__: March 16, 2017

### New Features and Changes

*	 Added ability to delete an entity.  Read more about this [here](http://docs.gdc.cancer.gov/Data_Submission_Portal/Users_Guide/Data_Upload_UG/#deleting-submitted-entities) <!-- SUBP-144 -->
*  Added Project Reports in the projects list page.  Read more about this [here](http://docs.gdc.cancer.gov/Data_Submission_Portal/Users_Guide/Homepage/#reports). <!-- SUBP-281, DAT-286, DAT-287, DAT-289 -->
*  To avoid confusion, renamed "Status" to "State" in the Browse section <!-- SUBP-237 -->
*  Added tooltip over Hierarchy title when reviewing an entity <!-- SUBP-238 -->
*  Restrict the upload window to only supported data formats (JSON and TSV) <!-- SUBP-245 -->

### Bugs Fixed Since Last Release

*  File status column should not be displayed for any clinical or biospecimen entities but only for submittable data files. <!--SUBP-226-->
*	Diagnosis / Treatment detail: Submitter ID (of the child / parent) is missing in the Details -> Hierarchy view. <!--SUBP-227-->
*	In some situations tooltip entries remain on-screen. Workaround is to refresh the page. <!--SUBP-229-->
*	In Browse tab, "Submittable Data Files" filter, clicking on "Download All" currently returns case and clinical informations instead of returning file informations. Workaround is to download information from the file the details panel. <!--SUBP-230-->
*	In Dashboard, the donut chart for number of cases with submittable data files is always empty. A workaround is to visit the Browse, detailed case view section to see, case by case, if it has submittable data files.<!--SUBP-231 / SUBP-156-->
*	In Transactions tab, after clicking on Commit or Discard, status is not automatically refreshed. <!--SUBP-232-->
*  Added the API version in the Data Submission Portal footer on the project list page. <!--SUBP-235-->
*  Inconsistent behavior when clicking on a Transaction ID on the Dashboard. <!--SUBP-239-->
*  Empty transactions created when submitting files in an incorrect format. <!--SUBP-242-->
*  JSON file downloaded from the Data Submission Portal cannot be used to resubmit data. <!--SUBP-243-->
*  "Submitted data files" donut chart and "Download Manifest" button do not get refreshed after committing a transaction. <!--SUBP-248-->
*  Release information on the Dashboard creates confusion. <!--SUBP-252-->
*  Missing Boolean fields from details panel. <!--SUBP-254-->
*  No message displayed if no results are found via the top menu search. <!--SUBP-256-->
*  "Invalid Date" on IE11 and Firefox ESR 45.x <!--SUBP-257-->
*  Download option truncated in the details panel. <!--SUBP-258-->
*  Download All from the browse section returns too many records. <!--SUBP-266-->

### Known Issues and Workarounds

*  When creating entities in the Submission Portal, occasionally an extra transaction will appear with status error. This does not seem to impact that actual transaction, which is recorded as occurring successfully.
<!--API-219-->
Release details are maintained in the [GDC Data Submission Portal Change Log](https://github.com/NCI-GDC/submission-ui/blob/master/CHANGELOG.md).

## Release 1.3.0

* __GDC Product__: GDC Data Submission Portal
* __Release Date__: October 31, 2016

### New Features and Changes

Not Applicable

### Bugs Fixed Since Last Release

*	Adding a call to the backend to ensure a refreshed token is being downloaded when user clicks on "Download Token". <!-- SUBP-170 -->
*	Fixed an issue with some buttons not working in Firefox <!-- SUBP-60 -->
*	Disabled Download Clinical button if the file has no clinical data <!-- SUBP-173 -->
*	Disabled Download Manifest button while the page is loading  <!-- SUBP-206 -->
* 	Fixed an issue with some dropdowns being cut off <!-- SUBP-212 / SUBP-214 -->

### Known Issues and Workarounds

*	Project submission and release is currently disabled. <!--SUBP-201-->
*	Diagnosis / Treatment detail: Submitter ID (of the child / parent) is missing in the Details -> Hierarchy view. <!--SUBP-227-->
* File status column should not be displayed for any clinical or biospecimen entities but only for submittable data files. <!--SUBP-226-->
*	In some situations tooltip entries remain on-screen. Workaround is to refresh the page. <!--SUBP-229-->
*	In Browse tab, "Submittable Data Files" filter, clicking on "Download All" currently returns case and clinical informations instead of returning file informations. Workaround is to download information from the file the details panel. <!--SUBP-230-->
*	In Dashboard, the donut chart for number of cases with submittable data files is always empty. A workaround is to visit the Browse, detailed case view section to see, case by case, if it has submittable data files.<!--SUBP-231 / SUBP-156-->
*	In Transactions tab, after clicking on Commit or Discard, status is not automatically refreshed. Workaround is to refresh the page after clicking on Commit or Discard. This does not affect the transaction section of the project dashboard. <!--SUBP-232-->
*  Reports are currently not available in the Data Submission Portal and will be added back in an upcoming version:
    *   Data Validation Report: The rows in the report are sometimes duplicated and #Files in error are not showing up in the report. The user should go to Project > Browse > Submitted Files to see the files in error and the error type<!--PGDC-1025 and PGDC-997-->.
    *   The Scientific Pre-alignment QC Report is not available.

Release details are maintained in the [GDC Data Submission Portal Change Log](https://github.com/NCI-GDC/submission-ui/blob/master/CHANGELOG.md).

## Release 1.2.2

* __GDC Product__: GDC Data Submission Portal
* __Release Date__: September 23, 2016

### New Features and Changes

This version contains major improvements to the GDC Data Submission Portal in both usability, performance and reliability.

Some known issues and workarounds listed in previous release notes have been made redundant due to this refactoring effort, thus are not listed anymore.

Please refer to the GDC Data Submission Portal User's Guide for more details about the features.

* Submission-related actions have been made Asynchronous.
* Fully revamped the dashboard layout and features to clarify the submission process and give easier access to key features.
* Created a transactions list page with options to take actions on transactions (in particular committing an upload)
* Improved performance of the Browse tab.
* Added GDC Apps to the header section.

### Bugs Fixed Since Last Release

* Data submitted to the project can be downloaded from each project page by clicking on "PROJECT DATA" from the project page.<!--PGDC-774-->
* When uploading multiple files at once, validation will fail if a child entity is listed before its parent. <!--PGDC-861-->
* In Browse > Case > Details, Experimental Data (renamed to Submittable Data Files) are not listed in "Related Entities" section. <!--PGDC-850-->
* In the upload report, the number of affected cases is incorrect (show 0) when entities are created. <!--PGDC-838-->


### Known Issues and Workarounds

*	Project submission and release is currently disabled. <!--SUBP-201-->
*	If case has no clinical data, the "Download Clinical" button is not disabled. The downloaded TSV will not contain Clinical Data. <!--SUBP-173-->
*	Download Manifest button is available while page is loading or when no files are in "registered" state. Clicking on Download will return a file with an error message. <!--SUBP-206-->
*	In Internet Explorer, GDC APPs and File dropdown are incorrectly aligned, making some elements only partially visible. <!--SUBP-212--><!--SUBP-214-->
*	Diagnosis / Treatment detail: Submitter ID (of the child / parent) is missing in the Details -> Hierarchy view. <!--SUBP-227-->
* File status column should not be displayed for any clinical or biospecimen entities but only for submittable data files. <!--SUBP-226-->
*	In some situations tooltip entries remain on-screen. Workaround is to refresh the page. <!--SUBP-229-->
*	In Browse tab, "Submittable Data Files" filter, clicking on "Download All" currently returns case and clinical informations instead of returning file informations. Workaround is to download information from the file the details panel. <!--SUBP-230-->
*	In Dashboard, the donut chart for number of cases with submittable data files is always empty. A workaround is to visit the Browse, detailed case view section to see, case by case, if it has submittable data files.<!--SUBP-231-->
*	In Transactions tab, after clicking on Commit or Discard, status is not automatically refreshed. Workaround is to refresh the page after clicking on Commit or Discard. This does not affect the transaction section of the project dashboard. <!--SUBP-232-->
*  Reports are currently not available in the Data Submission Portal and will be added back in an upcoming version:
    *   Data Validation Report: The rows in the report are sometimes duplicated and #Files in error are not showing up in the report. The user should go to Project > Browse > Submitted Files to see the files in error and the error type<!--PGDC-1025 and PGDC-997-->.
    *   The Scientific Pre-alignment QC Report is not available.

Release details are maintained in the [GDC Data Submission Portal Change Log](https://github.com/NCI-GDC/submission-ui/blob/master/CHANGELOG.md).

## Release 1.1.0

* __GDC Product__: GDC Data Submission Portal
* __Release Date__: May 20th, 2016

### New Features and Changes

* Updated login text

### Bugs Fixed Since Last Release

* Improved 508 compliance of the landing page

### Known Issues and Workarounds

*	Some actions on the data submission portal are resource intensive and might result in timeout or errors. The team is currently working on a major update to the data submission portal expected to be available towards the end of the summer. This version will address performance issues and improve overall user experience.

Release details are maintained in the [GDC Data Submission Portal Change Log](https://github.com/NCI-GDC/submission-ui/blob/master/CHANGELOG.md).

## Release 0.3.24.1

* __GDC Product__: GDC Data Submission Portal
* __Release Date__: February 26, 2016

### New Features and Changes

* Fully revamped dashboard, improved widgets with more details.
* Improved the release process and provided guidance.
* Removed dictionary viewer (moved to Documentation site).
* Added reports.

### Bugs Fixed Since Last Release

* Issues with transactions list page have been addressed.
* Search feature has been improved.
* Release processes (submit and release) have been finalized and implemented.

### Known Issues and Workarounds

* In Browse > Submitted Files, "DOWNLOAD ALL" downloads cases and clinical data instead of submitted files. A workaround is to initiate the download from the Dashboard<!--PGDC-970-->.
* In Browse > Case > Details, experimental data are not listed in "Related Entities" section<!--PGDC-850-->.
* When uploading multiple files at once, validation will fail if a child entity is listed before its parent. A workaround is to ensure parent entities are listed before their children in the Upload wizard<!--PGDC-861-->.
* Submission status column is inconsistent between Submitted Files and Read Group. Submission Status says "Validated" when it should say "Submitted" when users submit data to the GDC<!--PGDC-959-->.
* In the upload report, the number of affected cases is incorrect (show 0) when entities are created<!--PGDC-838-->.
* In Dashboard > Release, Download of submitted data returns data from the project workspace but not from snapshot (submitted files)<!--PGDC-774-->.
* In Browse > Read Groups, the data completeness property is incorrect when multiple files are in the bundle<!--PGDC-916-->.
* In Browse > Diagnosis/Treatment > Details, the hierarchy section is missing elements<!--PGDC-1023-->.
* Data Validation Report: The rows in the report are sometimes duplicated and #Files in error are not showing up in the report. The user should go to Project > Browse > Submitted Files to see the files in error and the error type<!--PGDC-1025 and PGDC-997-->.
* The Scientific Pre-alignment QC Report is not available.
* The Submission Portal is not meant to support XML file submission, users have to submit files through the GDC API.

Release details are maintained in the [GDC Data Submission Portal Change Log](https://github.com/NCI-GDC/submission-ui/blob/master/CHANGELOG.md).

## Release 0.3.21

* __GDC Product__: GDC Data Submission Portal
* __Release Date__: January 27, 2016

### New Features and Changes

*   New landing page.
*   Updated dashboard with new widgets, updated instructions and added a transactions list.
*   Improved user experience in the browse section of the portal (new layout, revamped detailed section, updated icons, etc.).

### Bugs Fixed Since Last Release

*   Some sections of the GDC Data Submission Portal are not compatible with Internet Explorer 10 and 11 (Note: support has been dropped for IE 10).
*   Case count coverage report is currently not available. Case overview report is available from the landing page.
*   If a page does not return data, the Submission Portal does not indicate that there is no data.
*   The projects dropdown show legacy projects. It should only show projects active for submission.

### Known Issues and Workarounds

*   Users may experience issues (e.g., unable to access a specific transaction ID) when navigating the transactions list page. The workaround is to refresh the page.
*   Search feature partially implemented, some items are not searchable.
*   Implementation of Submit and Release not finalized.
*   XML file submission is returning a bad request error. The workaround is to submit files through the GDC API.
*   Data Bundles - Lane Level Sequence: read group ID not unique, the message generated is not user friendly.

Release details are maintained in the [GDC Data Submission Portal Change Log](https://github.com/NCI-GDC/submission-ui/blob/master/CHANGELOG.md).

## Release 0.2.18.3

* __GDC Product__: GDC Data Submission Portal
* __Release Date__: November 30, 2015

### New Features and Changes

*   Dashboard to provide high-level details about one or more projects.
*   Submission wizard to guide user through a three stage submission:
    *   Upload file to the GDC.
    *   Validate uploaded files.
    *   Submit validated files to the GDC.
*   Support submission of JSON, TSV and XML files.
*   Tables and reports to identify key elements of the submission:
    *   List of Cases, Samples, Portions, Analytes, Aliquots, Lane Level Sequence.
    *   Cases missing clinical data.
    *   Cases missing samples.
    *   Aliquots missing data bundles.
    *   Submitted data bundle and associated QC metrics.
    *   Recent and all transactions.
    *   Case count coverage.
*   Manifest to facilitate submission of molecular data via the GDC Data Transfer Tool.
*   Per-project dictionary viewer.
*   Download authentication token.

### Bugs Fixed Since Last Release

*   Initial Release - Not Applicable

### Known Issues and Workarounds

*   Some sections of the GDC Data Submission Portal are not compatible with Internet Explorer 10 and 11.
*   Case count coverage report is currently not available.
*   Unable to use the search feature (implementation of this feature is not complete).
*   If the page does not return data, the Submission Portal does not mention that there is no data.
*   Case is not releasable when submitting a partial lane level seq data bundle.
    <u>Note</u>: This feature is currently being revised.
*   XML file submission is returning bad request error.
    A workaround is to submit files through the API.
*   Data Bundles - Lane Level Sequence: read group ID not unique, message generated is not user friendly.
*   The projects dropdown show legacy projects, it should only show projects active for submission.

Release details are maintained in the [GDC Data Submission Portal Change Log](https://github.com/NCI-GDC/submission-ui/blob/master/CHANGELOG.md).
