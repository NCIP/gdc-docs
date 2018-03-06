# Data Dictionary Release Notes

| Version | Date |
|---|---|
| v.1.11 | January 20, 2018 |
| v1.10.0 | August 22, 2017 |
| v1.7.1 | March 16, 2017 |
| v1.3.1 | September 7, 2016 |

---

## v.1.11

* __GDC Product__: GDC Data Dictionary
* __Release Date__: January 20, 2018


### New Features and Changes

* Added a link between the `sample` and `analyte` entities in the data model <!--TT-260--> <!--TT-259-->
* Added a link between the `sample` entity and other `sample` entities in the data model <!--TT-261-->
* Created `structural_variant_calling_workflow` entity <!--TT-270-->
* Created `structural_variation` entity <!--TT-207-->
* Removed `clinical_test` entity  <!--TT-232-->
* Removed `exon_expression` entity <!--TT-203-->
* Modified relationships of entities
    * Changed relationship between `submitted_genomic_profile` and `read_group` to one-to-many <!--TT-75-->
    * Changed relationship between `projects` and `masked_somatic_mutations` from one-to-one to many_to_one <!--TT-172-->
    * Changed relationship between `projects` and `aggregated_somatic_mutations` from one-to-one to many_to_one <!--TT-173-->
    * Changed relationship of `rna_expression_workflow` to downstream entities from one-to-one to one-to-many <!--TT-171-->
    * Changed relationship of `alignment_workflow` to downstream entities from many-to-one to many-to-many <!--TT-174-->
* Modified `project` entity
    - Added new state
        - `processed` <!--TT-266-->
    - Added new fields
        - `release_requested` <!--TT-267-->
        - `awg_review` <!--TT-267-->
        - `is_legacy` <!--TT-252-->
* Modified `clinical_supplement` entity
    - Added new `data_format` field <!--TT-271-->
        - `CDC JSON`
* Modified `case` entity
    - Enumerated `primary_site` field <!--TT-189-->
* Modified `diagnosis` entity
    - Enumerated fields
        - `primary_diagnosis` <!--TT-185-->
        - `tissue_or_organ_of_origin` <!--TT-190-->
        - `site_of_resection_or_biopsy` <!--TT-192-->
        - `tumor_grade` <!--TT-193-->
        - `morphology` <!--TT-186-->
        - `ajcc_clinical_stage` <!--TT-238-->
        - `ajcc_pathologic_stage` <!--TT-239-->
        - `laterality` <!--TT-242-->
        - `method_of_diagnosis` <!--TT-244-->
    - Changed type of fields
        - `year_of_diagnosis`: changed type to `int`<!--TT-121-->
        - `age_at_diagnosis`: changed type to `int` <!--TT-121-->
    - Added new permissible values to fields
        - `figo_stage`
            - `Stage IIC` <!--TT-241-->
        - `lymphatic_invasion_present`
            - `Not Reported` <!--TT-243-->
        - `perineural_invasion_present`
            - `Not Reported` <!--TT-245-->
        - `ann_arbor_clinical_stage`
            - `Not Reported` <!--TT-222-->  
            - `Unknown`<!--TT-240-->
        - `residual_disease`<!--TT-265-->
            - `Not Reported`
            - `Unknown`
        - `vascular_invasion_type`<!--TT-246-->
            - `Macro`
            - `Micro`
            - `No Vascular Invasion`
            - `Not Reported`
            - `Unknown`    
* Modified `exposure` entity
    - Enumerated fields
        - `alcohol_intensity` <!--TT-191-->
        - `alcohol_history` <!--TT-247-->
        - `cause_of_death` <!--TT-264-->
    - Removed `6` as a valid field from `tobacco_smoking_status` <!--TT-214-->
* Modified `demographic` entity
    - Enumerated field
        - `comorbidity` <!--TT-290-->
    - Added new permissible values to field
        - `cause_of_death` <!--TT-264-->
            - `Infection`
            - `Toxicity`
            - `Spinal Muscular Atrophy`
            - `End-stage Renal Disease`
* Modified `family history` entity
    - Enumerated fields
        - `relationship_primary_diagnosis` <!--TT-67-->
        - `relationship_type` <!--TT-67-->
* Modified `treatment` entity
    - Enumerated fields
        - `treatment_anatomic_site` <!--TT-221-->
        - `treatment_type` <!--TT-226-->
    - Changed type of fields
        - `days_to_treatment`: changed type to `int`<!--TT-213-->
        - `days_to_treatment_end`: changed type to `int`<!--TT-213-->
        - `days_to_treatment_start`: changed type to `int`<!--TT-213-->
* Modified `follow_up` entity
    - Enumerated field
        - `comorbidity` <!--TT-225-->
    - Removed `None` as valid field from `comorbidity` <!--TT-290-->
* Modified `demographic` entity
    - Added new field
        - `age_at_index` <!--TT-181-->
* Modified read_group entity <!--TT-200-->
    - Added new fields
        - `fragment_minimum_length`
        - `fragment_maximum_length`
        - `fragment_mean_length`
        - `fragment_standard_deviation_length`
* Modified `slide image` entity <!--TT-262-->
    - Added new data_formats
        - `JPEG` and `TIFF`
    - Added new data_type
        - `Cell Culture`
    - Added new experimental strategy
        - `Cell Culture`
    - Added new property
        - `Magnification`
    - Added new field
        - `date_time`
* Modified `sample` entity
    - Added new permissible values to fields <!--TT-274-->
        - `method_of_sample_procurement`
            - `Autopsy`  
        - `sample_type` <!--TT-257-->
            - `2D Classical Conditionally Reprogrammed Cells`
            - `2D Modified Conditionally Reprogrammed Cells`
            - `3D Organoid`
            - `3D Air-Liquid Interface Organoid`
            - `3D Neurosphere`
            - `Adherent Cell Line`
            - `Liquid Suspension Cell Line`
        - `composition` <!--TT-258-->
            - `2D Classical Conditionally Reprogrammed Cells`
            - `2D Modified Conditionally Reprogrammed Cells`
            - `3D Organoid`
            - `3D Air-Liquid Interface Organoid`
            - `3D Neurosphere`
            - `Adherent Cell Line`
            - `Liquid Suspension Cell Line`
* Modified `genomic_profile_harmonization_workflow` entity
    - Added new permissible values to field
        - `workflow_type` <!--TT-117-->
            - `GENIE Simple Somatic Mutation`
            - `GENIE Copy Number Variation`
* Modified `somatic_mutation_calling_workflow` entity
    - Added new permissible values to field
        - `workflow_type`<!--TT-206-->
            - `Pindel`
* Modified `rna_expression_workflow` entity <!--TT-202--> <!--TT-208-->
    - Added new permissible values to field
        - `workflow_type`
            - `RSEM - Quantification`
            - `STAR - Counts`
            - `RNA-SeQC - Counts`
            - `RNA-SeQC - FPKM`
            - `Kallisto - Quantification`
            - `Kallisto - HDF5`
* Modified `somatic_aggregation_workflow` <!--TT-202--> <!--TT-208-->
    - Added new permissible values to field
        - `workflow_type`
            - `Pindel Variant Aggregation and Masking`
            - `GENIE Variant Aggregation and Masking`
* Modified `gene_expression entity` <!--TT-204-->
    - Added new data_types
        - `Isoform Expression Quantification`
        - `Exon Expression Quantification`
    - Added new data_format
        - `HDF5`


### Bugs Fixed Since Last Release
* Fixed issue when submitting `tobacco_smoking_status` via tsv <!--TT-205-->

### Known Issues and Workarounds



## Release with API v1.10.0

* __GDC Product__: GDC Data Dictionary
* __Release Date__: August 22, 2017


### New Features and Changes

* Created `follow_up` entity to support longitudinal clinical data <!--TT-65-->
* Deprecated `clinical_test` entity <!--DOC-67-->
* Modified acceptable values for Read Group properties <!--TT-9,TT-76-->
    - `library_selection` : "Hybrid Selection, Affinity Enrichment, Poly-T Enrichment, Random, rRNA Depletion, miRNA Size Fractionation, Targeted Sequencing"
    - `library_strategy` : "Targeted Sequencing"
* Modified Diagnosis entity <!--TT-64-->
    - Added field `iss_stage`
    - Added field `best_overall_response`
    - Added field `days_to_best_overall_response`
    - Added field `progression_free_survival`
    - Added field `progression_free_survival_event`
    - Added field `overall_survival`
    - Added field `days_to_diagnosis` <!--TT-91-->
* Modified Treatment entity <!--TT-66-->
    - Added field `regimen_or_line_of_therapy`
* Modified Demographic entity <!--TT-83-->
    - Added field `cause_of_death`
    - Added field `days_to_birth`
    - Added field `days_to_death`
    - Added field `vital_status`
* Modified Case entity <!--TT-84-->
    - Added field `days_to_lost_to_followup`
    - Added field `lost_to_followup`
    - Added field `index_date`
* Added new tumor code, tumor id, and sample types to Sample entity to support OCG <!--TT-85, TT-68-->
    - `tumor_code` : "Acute leukemia of Ambiguous Lineage (ALAL), Lymphoid Normal, Tumor Adjacent Normal - Post Neo Adjuvant Therapy"
    - `tumor_code_id` : "15, 17, 18"
* Created `somatic_mutation_index` entity <!--TT-92-->
* Updated caDSR CDE links in data dictionary <!--DAT-794-->
* Added new `sample_type` : `tumor` to sample entity <!--TT-77-->
* Made `classification_of_tumor` on diagnosis entity non-required <!--DAT-203-->
* Added support for FM-AD to Genomic Profile Harmonization Workflow entity <!--DAT-985-->
* Added `data_type` : `Gene Level Copy Number Scores` to Copy Number Segment entity <!--TT-94-->

### Bugs Fixed Since Last Release

None

### Known Issues and Workarounds

*  Portion `weight` property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->

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
