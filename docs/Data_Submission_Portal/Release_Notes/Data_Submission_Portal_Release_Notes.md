# Data Submission Portal Release Notes

## Release 1.2.0

* __GDC Product__: GDC Data Submission Portal
* __Release Date__: September 22nd, 2016

### New Features and Changes

This version contains major improvements to the GDC Data Submission Portal in both usability, performance and reliability.

Some known issues and workaround listed in previous release notes have been made redundant due to this refactoring effort, thus are not listed anymore.

Please refer to the GDC Data Submission Portal user guide for more details about the features.

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
*	If case has no clinical data, the "Download Clinical" button is not disabled, the downloaded TSV will not contain Clinical Data. <!--SUBP-173-->
*	Download Manifest button is still available while page is loading, clickling on Download will return a file with an error message. <!--SUBP-206-->
*	In Internet Explorer, GDC APPs and File dropdown are incorrectly aligned, making some elements only partially visible. <!--SUBP-212--><!--SUBP-214-->
*	Diagnosis / Treatment detail: Submitter ID (of the child / parent) is missing in the Details -> Hierarchy view.  <!--SUBP-227-->
* File status column should not be displayed for any clinical or biospecimen entities but only for submittable data files. <!--SUBP-226-->
*	In some situtations tooltip entries remain on-screen. Workaround is to refresh the page. <!--SUBP-229-->
*	In Browse tab, "Submittable Data Files" filter, clicking on "Download All" currently returns case and clinical informations instead of returning file informations. Workaround is to download information from the file the details panel. <!--SUBP-230-->
*	In Dashboard, the donut chart for number of cases with submittable data files is always empty. A workaround is to visit the Browse, detailed case view section to see, case by case, if it has submittable data files.<!--SUBP-231-->
*	In Transactions tab, after clicking on Commit or Discard, status is not automatically refreshed. Workaround is to refresh the page after clicking on Commit or Discard. This does not affect the transaction section of the project dashboard. <!--SUBP-232-->
*  Reports are currently not available in the Data Submission Portal and will be added back in an upcoming version:
    *   Data Validation Report: The rows in the report are sometimes duplicated and #Files in error are not showing up in the report. The user should go to Project > Browse > Submitted Files to see the files in error and the error type<!--PGDC-1025 and PGDC-997-->.
    *   The Scientific Pre-alignment QC Report is not available.

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
