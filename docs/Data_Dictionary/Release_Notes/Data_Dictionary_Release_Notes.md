# Data Dictionary Release Notes


| Version | Date |
|---|---|
| [v.1.14](Data_Dictionary_Release_Notes.md#v114) | September 27, 2018 |
| [v1.13](Data_Dictionary_Release_Notes.md#v113) | May 21, 2018 |
| [v1.12.1](Data_Dictionary_Release_Notes.md#v1121) | April 26, 2018 |
| [v1.12](Data_Dictionary_Release_Notes.md#v112) | April 23, 2018 |
| [v1.11](Data_Dictionary_Release_Notes.md#v111) | January 20, 2018 |
| [v1.10.0](Data_Dictionary_Release_Notes.md#release-with-api-v1100) | August 22, 2017 |
| [v1.7.1](Data_Dictionary_Release_Notes.md#release-with-api-v171) | March 16, 2017 |
| [v1.3.1](Data_Dictionary_Release_Notes.md#release-with-api-v131) | September 7, 2016 |

## v.1.14

* __GDC Product__: GDC Data Dictionary
* __Release Date__: September 27, 2018

### New Features and Changes

* Added description and type information to generic properties <!--TT-1347-->
* Modify the dictionary viewer to be case sensitive for deprecated fields <!--PRTL-2180-->
* Deleted `sample_level_maf` entity  <!--TT-666-->
* Added new field for all nodes: <!--TT-556-->
    - `previous_versions_downloadable`
* Modified `molecular_test` entity
    - Migrated data from `blood_test` to `laboratory_test` and `biospecimen_type` <!--TT-754-->
    - Added new fields: <!--DAT-1644-->
        - `laboratory_test`
        - `biospecimen_type`
    - Added new permissible values for `gene_symbol` fields <!--DAT-1553-->
        - `Not Applicable`
    - Deleted field `blood_test` <!--DAT-1639-->
    - Add new permissible values for `antigen` field <!--DAT-1662-->
    - Add new permissible values to `molecular_analysis_method` <!--DAT-1663-->
    - Add new permissible values for `variant_type` field <!--DAT-1664-->
    - Add new permissible values to `test_result` <!--DAT-1665-->
* Modified `case` entity
    - Modified permissible values on `index_date`
        - Added new value `Initial Genomic Sequencing` <!--TT-1461-->
    - Updated CDE codes for properties <!--DAT-1590-->
    - Add new permissible value to `index_date` <!--DAT-1669-->
        - `Sample Procurement` <!--TT-678-->
* Modified `aliquot` entity
    - Changed `exclusive` to `True` <!--DAT-1519-->
* Modified `read_group_qc` entity
    - Made following fields not `required` and added `Unknown`/`Not Reported` as permissible values: <!--DAT-1519--><!--DAT-1525-->
        - `adapter_content`
        - `basic_statistics`
        - `encoding`
        - `kmer_content`
        - `overrepresented_sequences`
        - `per_base_sequence_quality`
        - `per_tile_sequence_quality`
        - `per_sequence_quality_score`
        - `per_base_sequence_content`
        - `per_sequence_gc_content`
        - `per_base_n_content`
        - `percent_gc_content`
        - `sequence_length_distribution`
        - `sequence_duplication_levels`
        - `total_sequences`
* Modified `aligned_reads` entity
    - Added new permissible values to `experimental_strategy` field <!--DAT-1526-->
        - `Targeted Sequencing`
        - `Bisulfite-Seq`
        - `ChIP-Seq`
        - `ATAC-Seq`
* Modified `read_group` entity
    - Added new field `days_to_sequencing` <!--DAT-1533-->
    - Add new permissible values to `target_capture_kit` <!--DAT-667-->
        - `Custom Personalis ACEcp VAREPOP-APOLLO Panel v2`
        - `xGen Exome Research Panel v1.0`
        - `SureSelect Human All Exon v5 + UTR`
* Modified `treatment` entity
    - Updated CDE codes for properties <!--DAT-1596-->
    - Updated description for `treatment_or_therapy` <!--DAT-1548-->
    - Updated description for `therapeutic_agents` <!--DAT-1549-->
    - Added new field:
        - `initial_disease_status` <!--DAT-1647-->
    - Updated permissible values for `treatment_type` <!--TT-672-->
    - Added new permissible values for `treatment_outcome`  <!--TT-680-->
        - `No Measurable Disease`
        - `Persistent Disease`
    - Deleted values from `treatment_intent_type`  <!--TT-681-->
* Modified `annotated_somatic_mutation` entity
    - Added new permissible values to `experimental_strategy` field <!--DAT-1554-->
        - `RNA-Seq`
        - `miRNA-Seq`
        - `Bisulfite-Seq`
        - `ChIP-Seq`
        - `ATAC-Seq`
    - Add new permissible value to `data_format` <!--TT-663-->
        - `MAF`
    - Add new permissible value to `data_type` <!--TT-665-->
        - `Masked Annotated Somatic Mutation`
* Modified `masked_somatic_mutation` entity
    - Added new permissible values to `experimental_strategy` field <!--DAT-1555-->
        - `RNA-Seq`
        - `miRNA-Seq`
        - `Bisulfite-Seq`
        - `ChIP-Seq`
        - `ATAC-Seq`
    - Removed `Validation` as permissible from `experimental_strategy` <!--TT-532-->
* Modified `simple_germlne_variation` entity
    - Added new permissible values to `experimental_strategy` field <!--DAT-1556-->
        - `RNA-Seq`
        - `miRNA-Seq`
        - `Bisulfite-Seq`
        - `ChIP-Seq`
        - `ATAC-Seq`
* Modified `simple_somatic_mutation` entity
    - Added new permissible values to `experimental_strategy` field <!--DAT-1558-->
        - `RNA-Seq`
        - `miRNA-Seq`
        - `Bisulfite-Seq`
        - `ChIP-Seq`
        - `ATAC-Seq`    
* Modified `submitted_genomic_profile` entity
    - Added new permissible values to `experimental_strategy` field <!--DAT-1558-->
        - `WXS`
        - `Low Pass WGS`
        - `RNA-Seq`
        - `miRNA-Seq`
        - `Bisulfite-Seq`
        - `ChIP-Seq`
        - `ATAC-Seq`
* Modified `diagnosis` entity
    - Updated CDE codes for properties <!--DAT-1592-->
    - Corrected the reference to `ubiquitous_properties` <!--DAT-1634-->
    - Deleted permissible values from `vascular_invasion_present` field <!--DAT-1560-->
        - `Extramural`
        - `Intramural`
    - Added permissible value to `ajcc_clinical_stage` <!--DAT-1571-->
        - `Stage IIIC1`
    - Added permissible values to `method_of_diagnosis` <!--DAT-1567-->
        - `Imaging`
        - `Physical Exam`
        - `Pathologic Review`
    - Added permissible value to `vascular_invasion_type` <!--DAT-1579-->
        - `Extramural`
        - `Intramural`  
    - Added permissible value to `metastasis_at_diagnosis_site`  <!--TT-673-->
    - Updated deprecated values list <!--TT-682-->
    - Deleted field `ldh_level` <!--TT-686-->
* Modified `follow_up` entity
    - Updated CDE codes for properties <!--DAT-1594-->
    - Added permissible value to `progression_or_recurrence_type` <!--DAT-1581-->
        - `Biochemical`
    - Added permissible values to `hpv_positive_type` <!--DAT-1583-->
        - `63`
        - `70`
    - Added permissible value to `disease_response` <!--DAT-1584-->
        - `BED-Biochemical Evidence of Disease`
        - `PDM-Persistent Distant Metastasis`
        - `PLD-Persistent Locoregional Disease`
        - `TF-Tumor Free`
        - `WT-With Tumor`
    - Added permissible value to `comorbidity` <!--DAT-1585-->
        - `Hepatitis A Infection`
        - `Herpes`      
    - Added permissible values to `risk_factor` field  <!--TT-674-->
* Modified `demographic` entity
    - Updated CDE codes for properties <!--DAT-1591-->
    - Modified permissible values for `cause_of_death` field   <!--TT-676-->
        - `Cardiovascular Disorder, NOS`
        - `Renal Disorder, NOS`
        - `Surgical Complications`
    - Updated `days_to` descriptions  <!--TT-683-->
* Modified `family_history` entity
    - Updated CDE codes for properties <!--DAT-1635-->
    - Added new fields: <!--DAT-1635-->
        - `relationship_primary_diagnosis`
        - `relationship_gender`
    - Updated CDE for field `relative_with_cancer_history` <!--TT-254-->
* Modified `submitted_unaligned_reads` entity
    - Make field `read_pair_number` required  <!--DAT-1648-->
    - Added permissible values to `read_pair_number` <!--DAT-1651-->
        - ` Not Applicable`
* Modified `aligned_reads` entity <!--DAT-1668-->
    - Added new fields:
        - `average_base_quality`
        - `average_insert_size`
        - `average_read_length`
        - `mean_coverage`
        - `pairs_on_diff_chr`
        - `proportion_base_mismatch`
        - `proportion_coverage_10X`
        - `proportion_coverage_30X`
        - `proportion_reads_duplicated`
        - `proportion_reads_mapped`
        - `proportion_targets_no_coverage`
        - `total_reads`
    - Changed nodes with `experimental_strategy` as `Validation` to `Targeted Sequencing` <!--TT-704-->
* Modified `sample` entity
    - Made `tissue_type` a required field <!--TT-657-->
    - Populated `tissue_type` with `Not Reported` for all nodes with no value  <!--TT-658-->
* Modified `alignment_workflow` entity
    - Added new field to `workflow_type` <!--TT-659-->
        - `STAR 2-Pass Chimeric`
* Modified `rna_expression_workflow` entity
    - Added new field to `workflow_type` <!--TT-660-->
        - `STAR - FPKM`
* Modified `gene_expression` entity
    - Add new permissible value to `data_type` <!--TT-663-->
        - `Splice Junction Quantification`
* Modified `aggregated_somatic_mutation` entity
    - Updated field `experimental_strategy` <!--TT-694-->
* Modified `somatic_mutation_calling_workflow` entity
    - Added `workflow_type` field <!--TT-661-->
* Modified `somatic_annotation_workflow` entity
    - Added `workflow_type` field <!--TT-662-->


### Bugs Fixed Since Last Release

* N/A

## v.1.13

* __GDC Product__: GDC Data Dictionary
* __Release Date__: May 21, 2018

### New Features and Changes

* Submitted Genomic Profile to Read Group relationship is updated to be many-to-many <!--DAT-1425-->
* Added requirement for submitter_id on additional nodes <!--TT-609-->

### Bugs Fixed Since Last Release

* None

## v.1.12.1

* __GDC Product__: GDC Data Dictionary
* __Release Date__: April 26, 2018

### New Features and Changes

* None

### Bugs Fixed Since Last Release

* Moved read_pair_number from Read Group to Submitted Unaligned Reads node <!--SV-1089-->

## v.1.12

* __GDC Product__: GDC Data Dictionary
* __Release Date__: April 23, 2018


### New Features and Changes

* Updates to data dictionary viewer:
    - Data dictionary viewer displays `true/false` instead of `boolean` <!--DOC-75-->
    - Updated `md5` description in data dictionary viewer to be more complete  <!--TT-400-->
* Creation of new entities:
    - `Copy Number Estimate` entity for copy number variation data <!--TT-439-->
    - `Copy Number Variation Workflow` entity for copy number variation pipeline metadata <!--TT-440-->
    - `Molecular Test` entity <!--TT-464-->
* Updates to all entities:
    - Updated all entities to include `batch_id` field for submission <!--TT-434-->
    - Updated all entities to include `downloadable` field <!--TT-501-->
    - Modified all file entities `file_size` field to be an integer <!--TT-457-->  <!--TT-403-->
    - Added support for creation of annotations on all entities in data model <!--DAT-1331-->
* Added new links / updated links between entities:
    - Created optional link between `follow_up` and `diagnosis` entities <!--TT-401-->
    - Created many-to-many link between `genomic_profile_harmonization_workflows` and `masked_somatic_mutations` <!--TT-530-->
    - Remove the restriction that prevents having alignment workflows to point to both `submitted_aligned_reads` and `submitted_unaligned_reads` <!--TT-422-->
* Updated the BCR XML endpoint to support submission of new entity relationships <!--DAT-1207-->
* Fixed description of `prior_malignancy` on `diagnosis` entity <!--TT-412-->
* Modified `sample` entity
    - Added new fields  <!--TT-399-->
        - `growth_rate`
        - `passage_count`
        - `catalog_reference`
        - `distributor_reference`
        - `distance_normal_to_tumor`
        - `biospecimen_laterality`
    - Added new permissible values to `composition` field  <!--TT-527-->
        - `Sorted Cells`
    - Added new permissible value to `sample_type` field <!--TT-509-->
        - `Next Generation Cancer Model`
    - Added new permissible values to `method_of_sample_procurement` field <!--TT-528-->
        - `Pancreatectomy`
        - `Whipple Procedure`
        - `Paracentesis`
    - Added new permissible values to `biospecimen_anatomic_site` field <!--TT-529-->
        - `Esophageal; Distal`
        - `Esophageal; Mid`
        - `Esophageal; Proximal`
        - `Hepatic Flexure`
        - `Rectosigmoid Junction`
    - Removed permissible values on `sample_type` <!--TT-510-->
        - `2D Classical Conditionally Reprogrammed Cells`
        - `2D Modified Conditionally Reprogrammed Cells`
        - `3D Organoid`
        - `3D Air-Liquid Interface Organoid`
        - `3D Neurosphere`
        - `Adherent Cell Line`
        - `Liquid Suspension Cell Line`
* Modified `aliquot` entity <!--TT-416-->
    - Added new fields
        - `no_matched_normal_wxs`
        - `no_matched_normal_wgs`
        - `no_matched_normal_targeted_sequencing`
        - `no_matched_normal_low_pass_wgs`
* Modified `read_group` entity
    - Updated descriptions on `read_group` entity  <!--TT-463-->
    - Added new field  <!--TT-415-->
        - `target_capture_kit`
    - Made `library_selection` a required field <!--TT-472-->
    - Removed duplicate `validators` field  <!--TT-484-->
    - Modified permissible values on `library_selection` field
        - Replaced `Affinity_Enrichment` with `Affinity Enrichment` <!--TT-447-->
        - Replaced `Poly-T_Enrichment` with `Poly-T Enrichment` <!--TT-449-->
        - Replaced `Hybrid_Selection` with `Hybrid Selection` <!--TT-451-->
        - Replaced `RNA_Depletion` with `rRNA Depletion` <!--TT-453-->
        - Replaced `Targeted Sequencing` with `Hybrid Selection` for FM-AD  <!--TT-505-->
    - Added new permissible value on `library_strand` field
        - `Not Applicable` <!--TT-402-->
    - Removed permissible values in `library_strategy` field <!--TT-442-->
        - `Amplicon`
        - `Validation`
        - `Other`
* Modified `diagnosis` entity
    - Added new fields  <!--TT-380-->
        - `ajcc_staging_system_edition`
        - `anaplasia_present`
        - `anaplasia_present_type`
        - `child_pugh_classification`
        - `cog_neuroblastoma_risk_group`
        - `cog_rhabdomyosarcoma_risk_group`
        - `enneking_msts_grade`
        - `enneking_msts_metastasis`
        - `enneking_msts_stage`
        - `enneking_msts_tumor_site`
        - `esophageal_columnar_dysplasia_degree`
        - `esophageal_columnar_metaplasia_present`
        - `first_symptom_prior_to_diagnosis`
        - `gastric_esophageal_junction_involvement`
        - `goblet_cells_columnar_mucosa_present`
        - `inpc_grade`
        - `inss_stage`
        - `irs_group`
        - `ishak_fibrosis_score`
        - `micropapillary_features`
        - `meduloblastoma_molecular_classification`
        - `metastasis_at_diagnosis`
        - `mitosis_karyorrhexis_index`
        - `peripancreatic_lymph_nodes_positive`
        - `peripancreatic_lymph_nodes_tested`
        - `synchronous_malignancy`
        - `supratentorial_localization`
        - `tumor_confined_to_organ_of_origin`
    - Added new permissible values to `vascular_invasion_present` field <!--TT-495-->
        - `Extramural`
        - `Intramural`
    - Updated `ann_arbor_clinical_stage` and `ann_arbor_pathologic_stage` fields  <!--TT-524-->
* Modified `aliquot` entity
    - Added new fields <!--TT-425-->
        - `selected_normal_wxs`
        - `selected_normal_wgs`
        - `selected_normal_targeted_sequencing`
        - `selected_normal_low_pass_wgs`
* Modified `project` entity
    - Added new fields
        - `request_submission`
        - `in_review` <!--TT-428-->
        - `submission_enabled` <!--TT-429-->
    - Removed duplicate `intended_release_date` field <!--TT-483-->
    - Made fields `not required` <!--TT-465-->
        - `disease_type`
        - `primary_site`
* Modified `case` entity
    - Updated descriptions  <!--TT-488-->
        - `disease_type`
        - `primary_site`
* Modified `demographic` entity
    - Added new fields <!--TT-379-->
        - `premature_at_birth`
        - `weeks_gestation_at_birth`
    - Made `submitter_id` a required field <!--TT-544-->
    - Changed type of `year_of_birth` from number to integer <!--TT-489-->
    - Changed type of `year_of_death` from number to integer <!--TT-490-->
* Modified `treatment` entity
    - Updated CDE code for `treatment_anatomic_site` <!--TT-492-->
    - Modified permissible values to `treatment_anatomic_site` field  <!--TT-493--> <!--TT-506-->
    - Made `submitter_id` a required field <!--TT-544-->
    - Added new permissible value to `treatment_intent_type` <!--TT-494-->
        - `Neoadjuvant`
    - Added new permissible value to `treatment_outcome` <!--TT-498-->
        - `Very Good Partial Response`
        - `Mixed Response`
        - `No Response`
    - Added new permissible value to `treatment_type` <!--TT-499-->
        - `Brachytherapy, High Dose`
        - `Brachytherapy, Low Dose`
        - `Radiation, 2D Conventional`
        - `Radiation, 3D Conformal`
        - `Radiation, Intensity-Modulated Radiotherapy`
        - `Radiation, Proton Beam`
        - `Radiation, Stereotactic Body`
        - `Radiation Therapy, NOS`
        - `Stereotactic Radiosurgery`
    - Removed field
        - `days_to_treatment`
    - Removed permissible values for `treatment_type` <!--TT-523-->
        - `Radiation`
        - `Radiation Therapy`
* Modified `analyte` entity
    - Removed duplicate `project` field <!--TT-486-->
* Modified `follow_up` entity
    - Added new permissible values to `comorbidity` field <!--TT-500-->
    - Added new fields <!--TT-365-->
        - `progression_or_recurrence_type`
        - `diabetes_treatment_type`
        - `reflux_treatment_type`
        - `barretts_esophagus_goblet_cells_present`
        - `karnofsky__performance_status`
        - `menopause_status`
        - `viral_hepatitis_serologies`
        - `reflux_treatment`
        - `pancreatitis_onset_year`
        - `comorbidity_method_of_diagnosis`
        - `risk_factor`
    - Removed fields <!--TT-517-->
        - `absolute_neutrophil`
        - `albumin`
        - `beta_2_microglobulin`
        - `bun`
        - `calcium`
        - `cea_level`
        - `colon_polyps_history`
        - `creatinine`
        - `crp`
        - `days_to_hiv_diagnosis`
        - `estrogen_receptor_percent_positive_ihc`
        - `estrogen_receptor_result_ihc`
        - `glucose`
        - `hemoglobin`
        - `her2_erbb2_percent_positive_ihc`
        - `her2_erbb2_result_fish`
        - `her2_erbb2_result_ihc`
        - `hiv_positive`
        - `hpv_status`
        - `iga`
        - `igg`
        - `igl_kappa`
        - `igl_lambda`
        - `igm`
        - `ldh_level`
        - `ldh_normal_range_upper`
        - `m_protein`
        - `microsatellite_instability_abnormal`
        - `platelet_count`
        - `progesterone_receptor_percent_positive_ihc`
        - `progesterone_receptor_result_ihc`
        - `total_protein`
        - `wbc`
* Modified `exposure` entity
    - Added new fields  <!--TT-383-->
        - `alcohol_days_per week`
        - `alcohol_drinks_per_day`
* Added `molecular_test` entity
    - Added new fields  <!--TT-384-->
        - `gene_symbol`
        - `second_gene_symbol`
        - `test_analyte_type`
        - `test_result`
        - `molecular_analysis_method`
        - `variant_type`
        - `molecular_consequence`
        - `chromosome`
        - `cytoband`
        - `exon`
        - `transcript`
        - `locus`
        - `dna_change`
        - `aa_change`
        - `rna_change`
        - `zygosity`
        - `histone_family`
        - `histone_variant`
        - `copy_number`
        - `antigen`
        - `test_value`
        - `test_units`
        - `specialized_molecular_test`
        - `ploidy`
        - `cell_count`
        - `loci_count`
        - `loci_abnormal_count`
        - `mismatch_repair_mutation`
        - `blood_test`
        - `blood_test_normal_range_upper`
        - `blood_test_normal_range_lower`
* Modified `biospecimen_supplement` entity
    - Added new permissible values to `data_format` field <!--TT-394-->
        - `TSV`
        - `BCR Biotab`
        - `BCR SSF XML`
        - `BCR PPS XML`
        - `FoundationOne XML`
        - `BCR Auxiliary XML`
    - Removed permissible values on `data_format` field <!--TT-397-->
        - `SSF`
        - `PPS`
* Modified `clinical_supplement` entity
    - Added new permissible values to `data_format` field <!--TT-395-->
        - `TSV`
        - `BCR OMF XML`
        - `BCR Biotab`
    - Removed field
        - `data_format`
* Modified `submitted_aligned_reads` entity
    - Added new permissible values to `experimental_strategy` field <!--TT-404-->
        - `Targeted Sequencing`
        - `Bisulfite-Seq`
        - `ChIP-Seq`
        - `ATAC-Seq`
    - Removed permissible values from `experimental_strategy` field <!--TT-404-->
        - `Validation`
        - `Total RNA-Seq`
* Modified `submitted_unaligned_reads` entity <!--TT-411-->
    - Added new field
        - `read_pair_number`
    - Added new permissible values to `experimental_strategy` field <!--TT-405-->
        - `Targeted Sequencing`
        - `Bisulfite-Seq`
        - `ChIP-Seq`
        - `ATAC-Seq`
* Modified `alignment_workflow` entity
    - Added new permissible values to `workflow_type` field  
        - `BWA with Mark Duplicates and BQSR` <!--TT-426-->
        - `BWA with BQSR` <!--TT-479-->
* Modified `copy_number_segment` entity
    - Removed permissible values for `data_type` field  <!--TT-441-->
        - `Gene Level Copy Number`
        - `Gene Level Copy Number Scores`
* Modified `submitted_genomic_profile` entity <!--TT-461-->
    - Added new permissible values for `data_type`
        - `Raw GCI Variant`
    - Added new permissible values for `data_category`
        - `Combined Nucleotide Variation`
    - Added new permissible values for `experimental_strategy`
        - `WGS`
    - Added new permissible values for `data_format`
        - `VCF`
* Modified `genomic_profile_harmonization_workflow` entity  <!--TT-466-->
    - Added new permissible value for `workflow_type`
        - `VCF LiftOver`
* Modified `simple_somatic_mutation` entity   <!--TT-467-->
    - Added new permissible values for `data_type`
        - `Raw GCI Variant`
    - Added new permissible values for `data_category`
        - `Combined Nucleotide Variation`
* Modified `masked_somatic_mutations` entity
    - Added new permissible value for `experimental_strategy` <!--TT-533-->
        - `Targeted Sequencing`
    - Removed permissible value for `experimental_strategy` <!--TT-532-->
        - `Validation`

### Bugs Fixed Since Last Release

* Fixed issue when `file_size` is specified as a float in submitted json file

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
