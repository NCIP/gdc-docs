# Data Dictionary Release Notes

## Release with API v1.10.0

* __GDC Product__: GDC Data Dictionary
* __Release Date__: August 22, 2017


### New Features and Changes

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

### Bugs Fixed Since Last Release

None

### Known Issues and Workarounds

*  Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->

## Release with API v1.7.1

* __GDC Product__: GDC Data Dictionary
* __Release Date__: March 16, 2017


### New Features and Changes

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

None

### Known Issues and Workarounds

None

## Release with API v1.3.1

* __GDC Product__: GDC Data Dictionary
* __Release Date__: September 7, 2016


### New Features and Changes

* Clinical Supplement entities can have `data_format` set to OMF.
* Biospecimen Supplement entities can have `data_format` set to SSF or PPS.
* Read group `instrument_model` can be set to "Illumina HiSeq 4000".
* Category of Slide entities in the GDC Data Model has changed from `data_bundle` to `biospecimen`.

### Bugs Fixed Since Last Release

None

### Known Issues and Workarounds

None
