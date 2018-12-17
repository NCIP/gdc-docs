# API Release Notes


| Version | Date |
|---|---|
| [v1.18.0](API_Release_Notes.md#v1160) | December 18, 2018 |
| [v1.17.0](API_Release_Notes.md#v1160) | November 7, 2018 |
| [v1.16.0](API_Release_Notes.md#v1160) | September 27, 2018 |
| [v1.15.0](API_Release_Notes.md#v1150) | August 23, 2018 |
| [v1.14.1](API_Release_Notes.md#v1141) | May 21, 2018 |
| [v1.14.0](API_Release_Notes.md#v1140) | April 23, 2018 |
| [v1.13.0](API_Release_Notes.md#v1130) | February 15, 2018 |
| [v1.12.0](API_Release_Notes.md#v1120) | December 21, 2017 |
| [v1.11.0](API_Release_Notes.md#v1110) | November 16, 2017 |
| [v1.10.0](API_Release_Notes.md#v1100) | August 22, 2017 |
| [v1.9.0](API_Release_Notes.md#v190) | June 29, 2017 |
| [v1.8.0](API_Release_Notes.md#v180) | May 9, 2017 |
| [v1.7.1](API_Release_Notes.md#v171) | March 16, 2017 |
| [v1.5.0](API_Release_Notes.md#v150) | October 31, 2016 |
| [v1.4.0](API_Release_Notes.md#v140) | September 23, 2016 |
| [v1.3.1](API_Release_Notes.md#v131) | September 7, 2016 |
| [v1.2.0](API_Release_Notes.md#v120) | August 9, 2016 |
| [v1.1.0](API_Release_Notes.md#v110) | May 25, 2016 |
| [v1.0.1](API_Release_Notes.md#v101) | May 16, 2016 |


## v1.18.0 <!--REQ-335-->

* __GDC Product__: Application Programming Interface (API)
* __Release Date__:  December 18, 2018

### New Features and Changes

* Update to auth for GDC Pre-Release Data Portal <!--TT-739-->

### Bugs Fixed Since Last Release

* None

### Known Issues and Workarounds

* Fields are not counted as missing if parent field is also missing.  This may occur with queries of nested fields in the Data Portal Advanced Search or an API query using a filter.  This behavior could impact results reported using search parameters of "IS MISSING" or "NOT MISSING". <!-- PGDC-2530 // https://github.com/NCI-GDC/gdcapi/pull/524  -->
* Certain very large API requests will time out.  It is recommended to break up very large requests into a series of smaller requests. <!-- PGDC-2411 -->



## v1.17.0

* __GDC Product__: Application Programming Interface (API)
* __Release Date__:  November 7, 2018

### New Features and Changes

* Created new index cnv_centric <!--TT743-->
* Created new index cnv_occurrence_centric <!--TT-744-->
* Created new REST API endpoints for CNV <!--TT-757-->
* Created mapping from aliquot to case for occurrence on cnv_centric <!--TT-762-->
* Created new graphql endpoints for CNV <!--TT-769-->
* Updated index case_centric to add cnv <!--TT-745-->
* Updated index gene_centric to add cnv <!--TT-747-->

### Bugs Fixed Since Last Release

* Fixed bug to prevent users from deleting files in state submitted or released <!--API-560-->

### Known Issues and Workarounds

* Fields are not counted as missing if parent field is also missing.  This may occur with queries of nested fields in the Data Portal Advanced Search or an API query using a filter.  This behavior could impact results reported using search parameters of "IS MISSING" or "NOT MISSING". <!-- PGDC-2530 // https://github.com/NCI-GDC/gdcapi/pull/524  -->
* Certain very large API requests will time out.  It is recommended to break up very large requests into a series of smaller requests. <!-- PGDC-2411 -->



## v1.16.0

* __GDC Product__: Application Programming Interface (API)
* __Release Date__: September 27, 2018

### New Features and Changes

* Tissue_type is now parsed from the BCR Biospecimen submission <!--TT-753-->


### Bugs Fixed Since Last Release

* Re-attach the /files path to the legacy blueprints in gdcapi <!--SV-1195-->
* Fix Internal Server Error when exporting JSON in Data Portal <!--SV-1209-->

### Known Issues and Workarounds

* Fields are not counted as missing if parent field is also missing.  This may occur with queries of nested fields in the Data Portal Advanced Search or an API query using a filter.  This behavior could impact results reported using search parameters of "IS MISSING" or "NOT MISSING". <!-- PGDC-2530 // https://github.com/NCI-GDC/gdcapi/pull/524  -->
* Certain very large API requests will time out.  It is recommended to break up very large requests into a series of smaller requests. <!-- PGDC-2411 -->




## v1.15.0

* __GDC Product__: Application Programming Interface (API)
* __Release Date__: August 23, 2018

### New Features and Changes

* Added support for file versioning <!--TT-482-->
* Added new API endpoints in support of versioning: history, files/versions, and files/versions/manifest  <!--PRTL-2015-->
* Updated the submit and release process <!--TT-249-->
* Created admin endpoint to allow for modifying file and project flags and states <!--API-368-->
* Populate version and release information for all existing files <!--TT-154-->

### Bugs Fixed Since Last Release

* None

### Known Issues and Workarounds

* Fields are not counted as missing if parent field is also missing.  This may occur with queries of nested fields in the Data Portal Advanced Search or an API query using a filter.  This behavior could impact results reported using search parameters of "IS MISSING" or "NOT MISSING". <!-- PGDC-2530 // https://github.com/NCI-GDC/gdcapi/pull/524  -->
* Certain very large API requests will time out.  It is recommended to break up very large requests into a series of smaller requests. <!-- PGDC-2411 -->


## v1.14.1

* __GDC Product__: Application Programming Interface (API)
* __Release Date__: May 21, 2018

### New Features and Changes

* None

### Bugs Fixed Since Last Release

* Update to XML submission to correctly import biospecimen relationnships <!--TT-578-->
* GraphQL submission queries return incorrect information <!--SV-1086-->

### Known Issues and Workarounds

* The error message for not including a token is incorrect for controlled access API requests (e.g. bam slicing).  It states that your token is invalid, but it should request that a token be supplied. <!--API-495-->
* Fields are not counted as missing if parent field is also missing.  This may occur with queries of nested fields in the Data Portal Advanced Search or an API query using a filter.  This behavior could impact results reported using search parameters of "IS MISSING" or "NOT MISSING". <!-- PGDC-2530 // https://github.com/NCI-GDC/gdcapi/pull/524  -->
* Certain very large API requests will time out.  It is recommended to break up very large requests into a series of smaller requests. <!-- PGDC-2411 -->


## v1.14.0

* __GDC Product__: Application Programming Interface (API)
* __Release Date__: April 23, 2018

### New Features and Changes

* Updated support of BCR XML submission for Next Generation Cancer Models <!--DAT-1207-->


### Bugs Fixed Since Last Release

* None

### Known Issues and Workarounds

* Fields are not counted as missing if parent field is also missing.  This may occur with queries of nested fields in the Data Portal Advanced Search or an API query using a filter.  This behavior could impact results reported using search parameters of "IS MISSING" or "NOT MISSING". <!-- PGDC-2530 // https://github.com/NCI-GDC/gdcapi/pull/524  -->
* Certain very large API requests will time out.  It is recommended to break up very large requests into a series of smaller requests. <!-- PGDC-2411 -->


## v1.13.0

* __GDC Product__: Application Programming Interface (API)
* __Release Date__: February 15, 2018

### New Features and Changes

* Reduced number of submission portal timeouts by optimizing graphql API requests <!--API-317-->

### Bugs Fixed Since Last Release

* Fixed bug preventing submission of clinical data via BCR XML <!--SV-985-->

### Known Issues and Workarounds

* Fields are not counted as missing if parent field is also missing.  This may occur with queries of nested fields in the Data Portal Advanced Search or an API query using a filter.  This behavior could impact results reported using search parameters of "IS MISSING" or "NOT MISSING". <!-- PGDC-2530 // https://github.com/NCI-GDC/gdcapi/pull/524  -->
* Certain very large API requests will time out.  It is recommended to break up very large requests into a series of smaller requests. <!-- PGDC-2411 -->

## v1.12.0

* __GDC Product__: Application Programming Interface (API)
* __Release Date__: December 21, 2017

### New Features and Changes

* None

### Bugs Fixed Since Last Release

* Updated submit project button on submission portal to work more quickly for projects with large numbers of nodes <!--API-344-->

### Known Issues and Workarounds

* Fields are not counted as missing if parent field is also missing.  This may occur with queries of nested fields in the Data Portal Advanced Search or an API query using a filter.  This behavior could impact results reported using search parameters of "IS MISSING" or "NOT MISSING". <!-- PGDC-2530 // https://github.com/NCI-GDC/gdcapi/pull/524  -->
* Certain very large API requests will time out.  It is recommended to break up very large requests into a series of smaller requests. <!-- PGDC-2411 -->



## v1.11.0

* __GDC Product__: Application Programming Interface (API)
* __Release Date__: November 16, 2017

### New Features and Changes

* Added support for histogram aggregation for numeric fields <!--API-355-->

### Bugs Fixed Since Last Release

* None

### Known Issues and Workarounds

* Fields are not counted as missing if parent field is also missing.  This may occur with queries of nested fields in the Data Portal Advanced Search or an API query using a filter.  This behavior could impact results reported using search parameters of "IS MISSING" or "NOT MISSING". <!-- PGDC-2530 // https://github.com/NCI-GDC/gdcapi/pull/524  -->
* Certain very large API requests will time out.  It is recommended to break up very large requests into a series of smaller requests. <!-- PGDC-2411 -->



## v1.10.0

* __GDC Product__: Application Programming Interface (API)
* __Release Date__: August 22, 2017

### New Features and Changes

* Longitudinal clinical data can now be submitted to the GDC
* Created new entities and reorganized the location of many properties in the GDC Dictionary.  A list of deprecated properties and their new locations can be found [here](../Release_Notes/DeprecatedFields_072017.txt).  Note that the deprecated properties still contain all information that they previously did.  In a later release we will map the data from the old location to the new.  Eventually, we will suspend submission to the old properties and will remove them entirely from the GDC.
* GDC Data Dictionary Changes Details
    * Created follow_up entity to support longitudinal clinical data <!--TT-65-->
    * Deprecated clinical test entity <!--DOC-67-->
    * Modified acceptable values for Read Group properties <!--TT-9,TT-76-->
    * Modified Diagnosis entity <!--TT-64-->  
    * Modified Treatment entity <!--TT-66-->
    * Modified Demographic entity <!--TT-83-->
    * Modified Case entity <!--TT-84-->
    * Added new tumor code, tumor id, and sample types to Sample entity to support OCG <!--TT-85, TT-68-->
    * Added property `days_to_diagnosis` to Diagnosis entity <!--TT-91-->
    * Created Somatic Mutation Index entity <!--TT-92-->
    * Updated CaDSR CDE links in data dictionary <!--DAT-794-->
    * Added new sample type `tumor` to sample entity <!--TT-77-->
    * Made classification_of_tumor on diagnosis entity non-required <!--DAT-203-->
    * Added support for FM-AD to Genomic Profile Harmonization Workflow entity <!--DAT-985-->
    * Added data type `Gene Level Copy Number Scores` to Copy Number Segment entity <!--TT-94-->



### Known Issues and Workarounds

* Fields are not counted as missing if parent field is also missing.  This may occur with queries of nested fields in the Data Portal Advanced Search or an API query using a filter.  This behavior could impact results reported using search parameters of "IS MISSING" or "NOT MISSING". <!-- PGDC-2530 // https://github.com/NCI-GDC/gdcapi/pull/524  -->
* Certain very large API requests will time out.  It is recommended to break up very large requests into a series of smaller requests. <!-- PGDC-2411 -->


## v1.9.0

* __GDC Product__: Application Programming Interface (API)
* __Release Date__: June 29, 2017

### New Features and Changes

* Updated ability to search data based on genes and mutations have been released. Discussion of the new endpoints and their use can be found in the Data Portal User Guide. These feature are limited to those projects for which the GDC contains open-access MAF files.
<!--API-112, API-114, API-116, API-118, API-119, APIT-179, API-180, API-181, API-181-->
* API request parameter `from` now returns first record when `from=0`.  Previously, `from=1` was the first record.  `from=0` is now the default.

### Known Issues and Workarounds

* Fields are not counted as missing if parent field is also missing.  This may occur with queries of nested fields in the Data Portal Advanced Search or an API query using a filter.  This behavior could impact results reported using search parameters of "IS MISSING" or "NOT MISSING". <!-- PGDC-2530 // https://github.com/NCI-GDC/gdcapi/pull/524  -->
* Certain very large API requests will time out.  It is recommended to break up very large requests into a series of smaller requests. <!-- PGDC-2411 -->


## v1.8.0

* __GDC Product__: Application Programming Interface (API)
* __Release Date__: May 9, 2017

### New Features and Changes

* Users can now control whether a set of files will be compressed or not when downloading.  For further details see the [API User Guide](../Users_Guide/Downloading_Files/#downloading-an-uncompressed-group-of-files). <!--API-175-->

### Known Issues and Workarounds

* Fields are not counted as missing if parent field is also missing.  This may occur with queries of nested fields in the Data Portal Advanced Search or an API query using a filter.  This behavior could impact results reported using search parameters of "IS MISSING" or "NOT MISSING". <!-- PGDC-2530 // https://github.com/NCI-GDC/gdcapi/pull/524  -->
* Certain very large API requests will time out.  It is recommended to break up very large requests into a series of smaller requests. <!-- PGDC-2411 -->


## v1.7.1

* __GDC Product__: Application Programming Interface (API)
* __Release Date__: March 16, 2017

### New Features and Changes

* Submission: Due to Data Dictionary updates new submission templates may be required for users submitting JSON and TSV formats
* Submission: Entities in submitted state (assigned when the project has been submitted) cannot be deleted.  <!-- API-129 -->
* Submission: When attempting to delete an entity that has child entities not specified in the request, an error message is generated that will include all of the child entities' UUIDs.  <!-- API-130 -->
* Submission: Entities associated with files uploaded to the GDC object store cannot be deleted until the associated file has been deleted. <!-- API-132 -->
* Re-enable Review, Submit, and Release functions for submission <!--API-73, API-138, API-159-->
* GDC Data Dictionary Changes
  * Added "submittable" property to all entities <!--DAT-215-->
  * Changed Read Group to category biospecimen <!--DAT-216-->
  * Added many new clinical properties available for submission <!--DAT-210, DAT-31, DAT-226, DAT-205-->
  * Added sample codes from Office of Cancer Genomics (OCG) to analyte and aliquot <!--DAT-170-->
  * Slides can now be attached to sample rather than just portion <!--DAT-205-->
  * `sample_type_id` is no longer required when submitting sample entities <!--DAT-233-->
  * `analyte_type_id` is no longer required when submitting aliquot and analyte entities<!--DAT-255-->
  * Clinical Test Entity is created for storing results of a variety of potential clinical tests related to the diagnosis - <!--DAT-223-->
  * Genomic Profiling Report entity created for storing particular derived sequencing results <!--DAT-229-->
  * Structural Variation entity created <!--DAT-229-->
  * Project entity includes new field "Intended Release Date" <!--API-143-->
  * Project entity includes new field "Releasable" <!--API-157-->


### Bugs Fixed Since Last Release

  * Fixed bug where boolean properties were not accepted with TSV submission <!--API-168-->

### Known Issues and Workarounds

* Fields are not counted as missing if parent field is also missing.  This may occur with queries of nested fields in the Data Portal Advanced Search or an API query using a filter.  This behavior could impact results reported using search parameters of "IS MISSING" or "NOT MISSING". <!-- PGDC-2530 // https://github.com/NCI-GDC/gdcapi/pull/524  -->
* Certain very large API requests will time out.  It is recommended to break up very large requests into a series of smaller requests. <!-- PGDC-2411 -->


## v1.5.0

* __GDC Product__: Application Programming Interface (API)
* __Release Date__: October 31, 2016

### New Features and Changes

* API responds with an error when the request specifies an unsupported combination of `filters` and `facets`. <!-- API-34 -->
* In TSV submissions, trailing and leading whitespace, including non-ASCII whitespace characters, are stripped from property names and values. <!-- API-63, API-29-->
* For released projects, any updates to previously submitted entities (i.e. `"state": "submitted"`) will be included in the following GDC data release. <!-- API-81 -->
* Performance improvements for manifest generation. <!-- API-58 -->

### Bugs Fixed Since Last Release

* Uploading certain unsupported metadata files caused the associated submission transactions to remain stuck in pending state. <!-- API-78 -->

### Known Issues and Workarounds

* Fields are not counted as missing if parent field is also missing.  This may occur with queries of nested fields in the Data Portal Advanced Search or an API query using a filter.  This behavior could impact results reported using search parameters of "IS MISSING" or "NOT MISSING". <!-- PGDC-2530 // https://github.com/NCI-GDC/gdcapi/pull/524  -->
* Certain very large API requests will time out.  It is recommended to break up very large requests into a series of smaller requests. <!-- PGDC-2411 -->




## v1.4.0

* __GDC Product__: Application Programming Interface (API)
* __Release Date__: September 23, 2016

### New Features and Changes

* Submission transaction log includes additional information to assist in tracking. <!-- API-40 -->
* Submission project state transitions are disabled temporarily while project release features are being improved. <!-- API-56 -->
* GDC data dictionary changes:
    * The `submittable` property was added to all entity types in the GDC data model. It indicates whether the entity type can be submitted by users. <!-- DM-16 -->
    * Category of Read Group entities in the GDC Data Model has changed from `data_bundle` to `biospecimen`.
    * Analyte entities support an expanded set of `analyte_type` values.

### Bugs Fixed Since Last Release

* None to report

### Known Issues and Workarounds

* API search & retrieval queries that do not include a `sort` parameter may return results in different order each time they are executed. This is a particular problem for paginated responses (i.e. responses to queries for which the number of results is greater than the `size` parameter). <!-- FEAT-120 -->
    * **Workaround:** Include a `sort` parameter in API search & retrieval queries.
* Fields are not counted as missing if parent field is also missing.  This may occur with queries of nested fields in the Data Portal Advanced Search or an API query using a filter.  This behavior could impact results reported using search parameters of "IS MISSING" or "NOT MISSING". <!-- PGDC-2530 // https://github.com/NCI-GDC/gdcapi/pull/524  -->
* Certain very large API requests will time out.  It is recommended to break up very large requests into a series of smaller requests. <!-- PGDC-2411 -->



## v1.3.1

* __GDC Product__: Application Programming Interface (API)
* __Release Date__: September 7, 2016

### New Features and Changes

* Successful `_dry_run` submission transactions can be committed to the GDC data model without having to re-upload metadata. The transactions can also be closed to prevent them from being committed in the future.
* Submission transactions can be submitted in asynchronous mode. In this mode, the GDC API will issue an immediate acknowledgement of the transaction, along with the `transaction_id`. The status of the transaction can be verified by the user at a later time by specifying the `transaction_id`. Users submitting large transactions may find this mode helpful.
* New submission transaction properties can be queried with GraphQL
* GDC Data Dictionary changes:
    * Clinical Supplement entities can have `data_format` set to OMF.
    * Biospecimen Supplement entities can have `data_format` set to SSF or PPS.
    * Read group `instrument_model` can be set to "Illumina HiSeq 4000".
    * Category of Slide entities in the GDC Data Model has changed from `data_bundle` to `biospecimen`.

### Bugs Fixed Since Last Release

* Incorrect BMI calculation in the import of BCR XML files.

### Known Issues and Workarounds

* API search & retrieval queries that do not include a `sort` parameter may return results in different order each time they are executed. This is a particular problem for paginated responses (i.e. responses to queries for which the number of results is greater than the `size` parameter). <!-- FEAT-120 -->
    * **Workaround:** Include a `sort` parameter in API search & retrieval queries.
* Fields are not counted as missing if parent field is also missing.  This may occur with queries of nested fields in the Data Portal Advanced Search or an API query using a filter.  This behavior could impact results reported using search parameters of "IS MISSING" or "NOT MISSING". <!-- PGDC-2530 // https://github.com/NCI-GDC/gdcapi/pull/524  -->
* Certain very large API requests will time out.  It is recommended to break up very large requests into a series of smaller requests. <!-- PGDC-2411 -->




## v1.2.0

* __GDC Product__: Application Programming Interface (API)
* __Release Date__: August 9, 2016

### New Features and Changes

* Tarballs generated by the `data` endpoint in response to multi-file data download requests now include a folder structure that puts each file in a folder whose name is the file's UUID. <!-- PGDC-2581 -->
* UUIDs in clinical XML files are no longer treated in a case-sensitive way by the `submission` endpoint. <!-- API-14 -->
* Improved performance of `submission` endpoint for transactions that include many cases. <!-- https://github.com/NCI-GDC/gdcapi/pull/495 -->
* Speed improvements for the `submission` endpoint.
* BCR XML is no longer validated against its XSD at submission.

### Bugs Fixed Since Last Release

* Fixed handling of `POST` requests to address problems with cart functionality in older versions of Firefox <!-- PGDC-2555 -->
* Files of category `related_files` can now be downloaded from the `data` endpoint. <!-- PGDC-2542 -->
* Allowed submission by projects in certain dbGaP registration states that were previously blocked. <!--API-3 -->

### Known Issues and Workarounds

* API search & retrieval queries that do not include a `sort` parameter may return results in different order each time they are executed. This is a particular problem for paginated responses (i.e. responses to queries for which the number of results is greater than the `size` parameter). <!-- FEAT-120 -->
    * **Workaround:** Include a `sort` parameter in API search & retrieval queries.
* Fields are not counted as missing if parent field is also missing.  This may occur with queries of nested fields in the Data Portal Advanced Search or an API query using a filter.  This behavior could impact results reported using search parameters of "IS MISSING" or "NOT MISSING". <!-- PGDC-2530 // https://github.com/NCI-GDC/gdcapi/pull/524  -->
* Certain very large API requests will time out.  It is recommended to break up very large requests into a series of smaller requests. <!-- PGDC-2411 -->




## v1.1.0

* __GDC Product__: Application Programming Interface (API)
* __Release Date__: May 25, 2016

### New Features and Changes

* BAM index files (.bai) are now automatically downloaded with parent BAM.

### Bugs Fixed Since Last Release

* None to report

### Bugs Fixed Since Last Release

* Sorting by file `submitter_id` no longer causes an internal server error
* BAM index files are now included with harmonized BAM files
* Certain very long API requests will time out.  It is recommended to break up into a series of smaller requests.

### Known Issues and Workarounds

* Fields are not counted as missing if parent field is also missing.  This may occur with queries of nested fields in the Data Portal Advanced Search or an API query using a filter.  This behavior could impact results reported using search parameters of "IS MISSING" or "NOT MISSING".
* Certain very large API requests will time out.  It is recommended to break up very large requests into a series of smaller requests.




## v1.0.1

* __GDC Product__: Application Programming Interface (API)
* __Release Date__: May 16, 2016

### New Features and Changes

* HTTP interface that uses JSON as the primary data exchange format
* Programmatic access to functionality provided by GDC Data and Submission portals, via `projects`, `cases`, `files`, `annotations`, `data`, `slicing`, `status`, and `submission` endpoints
* Programmatic access to GDC Legacy Archive via `legacy` endpoint
* Token-based authentication for secure access to controlled data and to submission functionality
* RESTful search that supports simple and complex queries via `filters`, `fields`, and `facets` parameters, and `project`, `files`, `cases`, and `annotations` endpoints.
* Search results can be sorted using `sort` parameter, paginated using `size` and `from` parameters, and output in JSON, TSV, and XML using `format` and `pretty` parameters.
* `_mapping` endpoint enables user discovery of fields available for data search and retrieval operations
* Support for downloading of individual files and of archives containing multiple files
* Generation of download and upload manifests for use with the GDC Data Transfer Tool
* BAM slicing functionality for downloading part(s) of a BAM file specified using chromosomal coordinates or HGNC gene names
* Transactional submission system that links individual data elements according to a graph-based GDC Data Model
* Two data entity identifiers: UUIDs, which are consistent across GDC, and Submitter IDs, for compatibility with submitters' tracking systems

### Bugs Fixed Since Last Release

* None to report


### Known Issues and Workarounds

* Use of non-ascii characters in token passed to Data Transfer Tool will produce incorrect error message "Internal server error: Auth service temporarily unavailable".
* Use of a decimal in an integer search field produces unexpected error.
* Certain very large API requests will time out.  It is recommended to break up very large requests into a series of smaller requests.
