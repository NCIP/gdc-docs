# Data Dictionary Release Notes

| Version | Date |
|---|---|
| [v.2.6.0](Data_Dictionary_Release_Notes.md#v260) | February 1, 2023 |
| [v.2.5.0](Data_Dictionary_Release_Notes.md#v250) | July 8, 2022 |
| [v.2.4.1](Data_Dictionary_Release_Notes.md#v241) | August 23, 2021 |
| [v.2.4.0](Data_Dictionary_Release_Notes.md#v240) | June 21, 2021 |
| [v.2.3.0](Data_Dictionary_Release_Notes.md#v230) | January 5, 2021 |
| [v.2.2.0](Data_Dictionary_Release_Notes.md#v220) | July 2, 2020 |
| [v.2.1.0](Data_Dictionary_Release_Notes.md#v210) | March 10, 2020 |
| [v.2.0.0](Data_Dictionary_Release_Notes.md#v200) | January 30, 2020 |
| [v.1.18.1](Data_Dictionary_Release_Notes.md#v1181) | November 6, 2019 |
| [v.1.18](Data_Dictionary_Release_Notes.md#v118) | July 31, 2019 |
| [v.1.17](Data_Dictionary_Release_Notes.md#v117) | June 5, 2019 |
| [v.1.16](Data_Dictionary_Release_Notes.md#v116) | April 17, 2019 |
| [v.1.15](Data_Dictionary_Release_Notes.md#v115) | December 18, 2018 |
| [v.1.14](Data_Dictionary_Release_Notes.md#v114) | September 27, 2018 |
| [v1.13](Data_Dictionary_Release_Notes.md#v113) | May 21, 2018 |
| [v1.12.1](Data_Dictionary_Release_Notes.md#v1121) | April 26, 2018 |
| [v1.12](Data_Dictionary_Release_Notes.md#v112) | April 23, 2018 |
| [v1.11](Data_Dictionary_Release_Notes.md#v111) | January 20, 2018 |
| [v1.10.0](Data_Dictionary_Release_Notes.md#release-with-api-v1100) | August 22, 2017 |
| [v1.7.1](Data_Dictionary_Release_Notes.md#release-with-api-v171) | March 16, 2017 |
| [v1.3.1](Data_Dictionary_Release_Notes.md#release-with-api-v131) | September 7, 2016 |

## v2.6.0

* __GDC Product__: GDC Data Dictionary
* __Release Date__: February 1, 2023

### New Features and Changes

* Altered `aliquot` Entity
	* New deprecated property: `analyte_type_id`
* Altered `analyte` Entity
	* New deprecated property: `analyte_type_id`
* Altered `annotation` Entity
	* Changes made to `links`
		* `copy_number_auxiliary_files` added to subgroup
	* New property: `copy_number_auxiliary_files`
* Altered `exposure` Entity
	* New deprecated property: `coal_dust_exposure`
	* New property: `alcohol_frequency`
	* New property: `chemical_exposure_type`
	* New property: `occupation_duration_years`
	* New property: `occupation_type`
	* Changes made to `exposure_type`
		* New permissable value: `Coal Dust`
	* Changes made to `tobacco_smoking_status`
		* New permissable value: `Current Reformed Smoker, Duration Not Specified`
		* New permissable value: `Current Reformed Smoker for < or = 15 yrs`
		* New permissable value: `Current Reformed Smoker for > 15 yrs`
		* New permissable value: `Current Smoker`
		* New permissable value: `Lifelong Non-Smoker`
		* New permissable value: `Smoker at Diagnosis`
		* New permissable value: `Smoking history not documented`
		* New deprecated value: `1`
		* New deprecated value: `2`
		* New deprecated value: `3`
		* New deprecated value: `4`
		* New deprecated value: `5`
		* New deprecated value: `6`
		* New deprecated value: `7`
* Altered `diagnosis` Entity
	* New property: `cancer_detection_method`
	* New property: `fab_morphology_code`
	* New property: `gleason_score`
	* New property: `pediatric_kidney_staging`
	* New property: `sites_of_involvement_count`
	* New property: `tumor_grade_category`
	* New property: `uicc_staging_system_edition`
	* New property: `weiss_assessment_findings`
	* Changes made to `inpc_grade`
		* New permissable value: `Undifferentiated or Poorly Differentiated`
	* Changes made to `sites_of_involvement`
		* Changes made to `items`
			* New permissable value: `Brain, Brain stem`
			* New permissable value: `Brain, Cerebellum, NOS`
			* New permissable value: `Brain, Cerebrum`
			* New permissable value: `Brain, Frontal lobe`
			* New permissable value: `Brain, Occipital lobe`
			* New permissable value: `Brain, Parietal`
			* New permissable value: `Brain, Temporal lobe`
			* New permissable value: `Brain, Ventricle, NOS`
			* New permissable value: `Central Lung`
			* New permissable value: `Gastrosplenic Ligament`
			* New permissable value: `Parametrium, Left`
			* New permissable value: `Parametrium, Right`
			* New permissable value: `Pelvic Lymph Node(s)`
			* New permissable value: `Peri-aortic Lymph Node(s)`
			* New permissable value: `Peripheral Lung`
			* New permissable value: `Thorax, NOS`
* Altered `family_history` Entity
	* Changes made to `relationship_type`
		* New permissable value: `First Degree Relative, NOS`
* Altered `follow_up` Entity
	* New deprecated property: `comorbidity`
	* New deprecated property: `risk_factor`
	* New property: `comorbidities`
	* New property: `days_to_first_event`
	* New property: `days_to_risk_factor`
	* New property: `evidence_of_progression_type`
	* New property: `first_event`
	* New property: `imaging_anatomic_site`
	* New property: `imaging_findings`
	* New property: `imaging_suv`
	* New property: `imaging_suv_max`
	* New property: `peritoneal_washing_results`
	* New property: `pregnancy_count`
	* New property: `risk_factors`
	* New property: `risk_factor_method_of_diagnosis`
	* New property: `timepoint_category`
	* New property: `year_of_follow_up`
	* Changes made to `evidence_of_recurrence_type`
		* New permissable value: `Histologic Confirmation`
	* Changes made to `history_of_tumor_type`
		* New permissable value: `Colorectal Cancer`
		* New permissable value: `Lower Grade Glioma`
	* Changes made to `pregnancy_outcome`
		* New permissable value: `Spontaneous Abortion`
	* Changes made to `progression_or_recurrence_type`
		* New permissable value: `Locoregional`
* Altered `gene_expression` Entity
	* Changes made to `experimental_strategy`
		* New permissable value: `m6A MeRIP-Seq`
* Altered `molecular_test` Entity
	* New property: `aneuploidy`
	* New property: `chromosomal_translocation`
	* New property: `chromosome_arm`
	* New property: `mutation_codon`
	* New property: `test_value_range`
	* New property: `timepoint_category`
	* Changes made to `antigen`
		* New permissable value: `Immunoglobulin, Cytoplasmic`
		* New permissable value: `Immunoglobulin, Surface`
	* Changes made to `laboratory_test`
		* New permissable value: `Blast Count`
		* New permissable value: `Metaphase Nucleus Count`
		* New permissable value: `Minimal Residual Disease`
		* New permissable value: `Monocytes`
	* Changes made to `test_result`
		* New permissable value: `Amplified`
		* New permissable value: `Elevated`
		* New permissable value: `Not Amplified`
	* Changes made to `variant_type`
		* New permissable value: `Fusion`
		* New permissable value: `Genetic Polymorphism`
		* New permissable value: `Internal Tandem Duplication`
		* New permissable value: `Mutation, NOS`
* Altered `pathology_detail` Entity
	* New property: `extracapsular_extension`
	* New property: `extranodal_extension`
	* New property: `extrascleral_extension`
	* New property: `lymph_node_dissection_method`
	* New property: `lymph_node_dissection_site`
	* New property: `lymph_nodes_removed`
	* New property: `percent_tumor_nuclei`
	* New property: `residual_tumor_measurement`
	* Changes made to `additional_pathology_findings`
		* New permissable value: `Bone marrow concordant histology`
		* New permissable value: `Bone marrow discordant histology`
	* Changes made to `dysplasia_type`
		* New permissable value: `Esophageal Mucosa Columnar Dysplasia`
	* Changes made to `lymph_node_involved_site`
		* New permissable value: `Aortic`
		* New permissable value: `Pelvis, NOS`
* Altered `read_group` Entity
	* Changes made to `target_capture_kit`
		* New property: `Custom SureSelect CGCI-BLGSP Panel - 7.8 Mb`
* Altered `sample` Entity
	* New required property: `preservation_method`
	* New required property: `specimen_type`
	* New required property: `tumor_descriptor`
	* New deprecated property: `composition`
	* New deprecated property: `sample_type`
	* New deprecated property: `sample_type_id`
	* New property: `specimen_type`
	* Changes made to `tumor_descriptor`
		* New permissable value: `New Primary`
* Altered `simple_germline_variation` Entity
	* Changes made to `experimental_strategy`
		* New permissable value: `Genotyping Array`
* Altered `somatic_copy_number_workflow` Entity
	* Changes made to `workflow_type`
		* New permissable value: `ABSOLUTE LiftOver`
		* New permissable value: `ASCAT3`
		* New permissable value: `GATK4 CNV`
* Altered `structural_variant_calling_workflow` Entity
	* Changes made to `workflow_type`
		* New permissable value: `SvABA`
* Altered `submitted_genotyping_array` Entity
	* Changed `aliquot` link from `one-to-one` to `many-to-one`
* Altered `treatment` Entity
	* New property: `residual_disease`
	* New property: `timepoint_category`
	* New property: `treatment_duration`
	* New property: `treatment_outcome_duration`
	* New property: `therapeutic_level_achieved`
	* New property: `therapeutic_target_level`
	* Changes made to `initial_disease_status`
		* New permissable value: `Persistent Disease`
		* New permissable value: `Refractory Disease`
	* Changes made to `protocol_identifier`
		* New permissable value: `901`
		* New permissable value: `911`
		* New permissable value: `914`
		* New permissable value: `925`
		* New permissable value: `935`
		* New permissable value: `937`
		* New permissable value: `3881`
		* New permissable value: `3891`
		* New permissable value: `4941`
		* New permissable value: `8605`
		* New permissable value: `9047`
		* New permissable value: `9082`
		* New permissable value: `9340`
		* New permissable value: `9341`
		* New permissable value: `9342`
		* New permissable value: `9343`
		* New permissable value: `9464`
		* New permissable value: `9640`
		* New permissable value: `9906`
		* New permissable value: `321P2`
		* New permissable value: `321P3`
		* New permissable value: `323P`
		* New permissable value: `A3961`
		* New permissable value: `A3973`
		* New permissable value: `AADM01P1`
		* New permissable value: `AALL0031`
		* New permissable value: `AALL0232`
		* New permissable value: `AALL0331`
		* New permissable value: `AALL03B1`
		* New permissable value: `AALL0434`
		* New permissable value: `AALL07P4`
		* New permissable value: `AALL08P1`
		* New permissable value: `AALL1131`
		* New permissable value: `AAML00P2`
		* New permissable value: `AAML03P1`
		* New permissable value: `AAML0531`
		* New permissable value: `AAML0631`
		* New permissable value: `AAML1031`
		* New permissable value: `AB9804`
		* New permissable value: `ACCL0331`
		* New permissable value: `ACCL0431`
		* New permissable value: `ACCL05C1`
		* New permissable value: `ACCL0934`
		* New permissable value: `ACCL1031`
		* New permissable value: `ADVL0018`
		* New permissable value: `ADVL0212`
		* New permissable value: `ADVL0214`
		* New permissable value: `ADVL0215`
		* New permissable value: `ADVL0421`
		* New permissable value: `ADVL0524`
		* New permissable value: `ADVL0525`
		* New permissable value: `ADVL06B1`
		* New permissable value: `ADVL0714`
		* New permissable value: `ADVL0812`
		* New permissable value: `ADVL0813`
		* New permissable value: `ADVL0821`
		* New permissable value: `ADVL0911`
		* New permissable value: `ADVL0912`
		* New permissable value: `ADVL0918`
		* New permissable value: `ADVL0921`
		* New permissable value: `ADVL1011`
		* New permissable value: `ADVL1111`
		* New permissable value: `ADVL1112`
		* New permissable value: `ADVL1115`
		* New permissable value: `ADVL1213`
		* New permissable value: `ADVL1412`
		* New permissable value: `AEPI07N1`
		* New permissable value: `ALTE03N1`
		* New permissable value: `ALTE05N1`
		* New permissable value: `ANBL0032`
		* New permissable value: `ANBL00B1`
		* New permissable value: `ANBL00P1`
		* New permissable value: `ANBL00P3`
		* New permissable value: `ANBL02P1`
		* New permissable value: `ANBL0321`
		* New permissable value: `ANBL0322`
		* New permissable value: `ANBL0421`
		* New permissable value: `ANBL0531`
		* New permissable value: `ANBL0532`
		* New permissable value: `ANBL0621`
		* New permissable value: `ANBL0931`
		* New permissable value: `ANBL09P1`
		* New permissable value: `ANBL1021`
		* New permissable value: `ANBL1221`
		* New permissable value: `ANUR1131`
		* New permissable value: `AOST0121`
		* New permissable value: `AOST0331`
		* New permissable value: `AOST06B1`
		* New permissable value: `AOST06P1`
		* New permissable value: `AREN03B2`
		* New permissable value: `B003`
		* New permissable value: `B903`
		* New permissable value: `B947`
		* New permissable value: `B954`
		* New permissable value: `B973`
		* New permissable value: `BCM`
		* New permissable value: `CCG2961`
		* New permissable value: `D9902`
		* New permissable value: `E04`
		* New permissable value: `E18`
		* New permissable value: `GBCTTO/99`
		* New permissable value: `GLATO 2006`
		* New permissable value: `I03`
		* New permissable value: `IHRT`
		* New permissable value: `INT-0133`
		* New permissable value: `N891`
		* New permissable value: `NWTS-4`
		* New permissable value: `NWTS-5`
		* New permissable value: `OSTEO 2006`
		* New permissable value: `Not Reported`
		* New permissable value: `P9462`
		* New permissable value: `P9641`
		* New permissable value: `P9754`
		* New permissable value: `P9851`
		* New permissable value: `P9906`
		* New permissable value: `P9963`
		* New permissable value: `R9702`
		* New permissable value: `S31`
		* New permissable value: `S921`
		* New permissable value: `STB`
	* Changes made to `treatment_anatomic_site`
		* New permissable value: `Locoregional Site`
	* Changes made to `treatment_intent_type`
		* New permissable value: `Consolidation Therapy`
		* New permissable value: `First-Line Therapy`
		* New permissable value: `Induction`
		* New permissable value: `Radiation Boost`
	* Changes made to `treatment_outcome`
		* New permissable value: `Normalization of Tumor Markers`
	* Changes made to `treatment_type`
		* New permissable value: `Biopsy, Excisional`
		* New permissable value: `Biopsy, Incisional`
		* New permissable value: `Distal Pancreatectomy`
		* New permissable value: `Surgery, NOS`
		* New permissable value: `Surgery, Minimally Invasive`
		* New permissable value: `Surgery, Open`
		* New permissable value: `Total Pancreatectomy`
		* New permissable value: `Whipple`
		* New deprecated value: `Cryoablation`
		* New deprecated value: `Surgery`
	* Changes made to `therapeutic_agents`
		* New permissable value: `TLR9 Agonist SD-101`
* New Entity: `copy_number_auxiliary_file`

### Known Issues and Workarounds

* The [GDC Data Dictionary Viewer](https://docs.gdc.cancer.gov/Data_Dictionary/viewer/) on the [GDC Documentation Site](https://docs.gdc.cancer.gov) does not currently display permissible values for array type properties. This does not impact submission of permissible values for these properties. Data submitters that would like to submit data for these properties can contact the GDC Helpdesk (support@nci-gdc.datacommons.io) for a list of permissible values for the affected properties. This will be addressed in a future release. 

## v2.5.0

* __GDC Product__: GDC Data Dictionary
* __Release Date__: July 8, 2022

### New Features and Changes

* Altered `diagnosis` Entity
	* New property: `double_expressor_lymphoma`
	* New property: `double_hit_lymphoma`
	* New property: `max_tumor_bulk_site`
	* New property: `uicc_clinical_m`
	* New property: `uicc_clinical_n`
	* New property: `uicc_clinical_stage`
	* New property: `uicc_clinical_t`
	* New property: `uicc_pathologic_m`
	* New property: `uicc_pathologic_n`
	* New property: `uicc_pathologic_stage`
	* New property: `uicc_pathologic_t`
	* Changes made to `ajcc_pathologic_n`
		* New permissible value: `N2mi`
	* Changes made to `ajcc_pathologic_t`
		* New permissible value: `T1c2`
	* Changes made to `classification_of_tumor`
		* New permissible value: `Prior primary`
		* New permissible value: `Synchronous primary`
	* Changes made to `metastasis_at_diagnosis_site`
		* New permissible value: `Spleen`
		* New permissible value: `Stomach`
	* Changes made to `morphology`
		* New permissible value: `8211/6`
		* New permissible value: `8460/6`
		* New permissible value: `8520/6`
		* New permissible value: `9715/3`
	* Changes made to `primary_diagnosis`
		* New permissible value: `Breast implant-associated anaplastic large cell lymphoma`
		* New permissible value: `Indolent T-cell lymphoproliferative disorder of gastrointestinal tract`
	* Changes made to `sites_of_involvement`
		* New permissible value: `Abdomen`
		* New permissible value: `Adrenal gland, NOS`
		* New permissible value: `Adnexa`
		* New permissible value: `Anus, NOS`
		* New permissible value: `Appendix`
		* New permissible value: `Bladder, NOS`
		* New permissible value: `Bone, NOS`
		* New permissible value: `Bone marrow`
		* New permissible value: `Brain, NOS`
		* New permissible value: `Breast, NOS`
		* New permissible value: `Bronchus`
		* New permissible value: `Buccal mucosa`
		* New permissible value: `Central nervous system`
		* New permissible value: `Cerebrospinal fluid`
		* New permissible value: `Chest wall`
		* New permissible value: `Colon, NOS`
		* New permissible value: `Diaphragm`
		* New permissible value: `Ear`
		* New permissible value: `Epididymis`
		* New permissible value: `Esophagus`
		* New permissible value: `Eye`
		* New permissible value: `Gallbladder`
		* New permissible value: `Heart`
		* New permissible value: `Kidney, NOS`
		* New permissible value: `Larynx, NOS`
		* New permissible value: `Leptomeninges`
		* New permissible value: `Liver`
		* New permissible value: `Lower lobe, lung`
		* New permissible value: `Lung, NOS`
		* New permissible value: `Lymph node, NOS`
		* New permissible value: `Mandible`
		* New permissible value: `Maxilla`
		* New permissible value: `Mediastinal soft tissue`
		* New permissible value: `Mesentery`
		* New permissible value: `Middle lobe, lung`
		* New permissible value: `Mouth`
		* New permissible value: `Nasopharynx`
		* New permissible value: `Nasal soft tissue`
		* New permissible value: `Neck`
		* New permissible value: `Ocular orbits`
		* New permissible value: `Oropharynx`
		* New permissible value: `Pancreas`
		* New permissible value: `Parotid gland`
		* New permissible value: `Pelvis`
		* New permissible value: `Pericardium`
		* New permissible value: `Perineum`
		* New permissible value: `Peri-orbital soft tissue`
		* New permissible value: `Peripheral blood`
		* New permissible value: `Peripheral nervous system`
		* New permissible value: `Peritoneum, NOS`
		* New permissible value: `Pleura`
		* New permissible value: `Prostate`
		* New permissible value: `Rectum`
		* New permissible value: `Retroperitoneum`
		* New permissible value: `Salivary gland`
		* New permissible value: `Sigmoid colon`
		* New permissible value: `Sinus`
		* New permissible value: `Skin`
		* New permissible value: `Small intestine`
		* New permissible value: `Soft tissue`
		* New permissible value: `Spleen`
		* New permissible value: `Stomach`
		* New permissible value: `Testes`
		* New permissible value: `Thyroid`
		* New permissible value: `Tongue`
		* New permissible value: `Tonsil`
		* New permissible value: `Trachea`
		* New permissible value: `Transplanted kidney`
		* New permissible value: `Transverse colon`
		* New permissible value: `Upper lobe, lung`
		* New permissible value: `Vagina`
		* New permissible value: `Vertebrae`
	* Changes made to `tumor_depth`
		* New property: `minimum`
* Altered `submitted_unaligned_reads` Entity
	* Changes made to `experimental_strategy`
		* New permissible value: `m6A MeRIP-Seq`
		* Removed permissible value: `m6A RNA Methylation`
* Altered `pathology_detail` Entity
	* New property: `tumor_level_prostate`
	* New property: `zone_of_origin_prostate`
	* Changes made to `additional_pathology_findings`
		* New permissible value: `Percent follicular component <= 10%`
		* New permissible value: `Percent follicular component > 10%`
* Altered `exposure` Entity
	* New property: `asbestos_exposure_type`
	* New property: `age_at_last_exposure`
	* New property: `exposure_source`
	* New property: `use_per_day`
	* Removed property: `marijuana_use_per_week`
	* Removed property: `smokeless_tobacco_quit_age`
	* Removed property: `tobacco_use_per_day`
	* Changes made to `alcohol_intensity`
		* New permissible value: `Social Drinker`
	* Added maximum to `exposure_duration_years`
	* Changes made to `exposure_type`
		* New permissible value: `Asbestos`
		* New permissible value: `Chemical`
		* New permissible value: `Radon`
		* New permissible value: `Respirable Crystalline Silica`
		* New permissible value: `Smokeless Tobacco`
* Altered `case` Entity
	* Added min/max to `days_to_consent`
* Altered `analyte` Entity
	* Changes made to `analyte_type`
		* New permissible value: `m6A Enriched RNA`
* Altered `read_group` Entity
	* Changes made to `library_strategy`
		* New permissible value: `m6A MeRIP-Seq`
		* Removed permissible value: `m6A RNA Methylation`
	* Changes made to `target_capture_kit`
		* New permissible value: `Custom SeqCap EZ BeatAML Panel - 12.5 Mb`
		* New permissible value: `Custom Twist MP2PRT-WT Panel - 52 Kb`
* Altered `aligned_reads` Entity
	* Changes made to `experimental_strategy`
		* New permissible value: `m6A MeRIP-Seq`
* Altered `follow_up` Entity
	* Changes made to `comorbidity`
		* New permissible value: `CNS Infection`
		* New permissible value: `EBV Lymphoproliferation`
		* New permissible value: `Hodgkin Lymphoma`
		* New permissible value: `Lymphamatoid Papulosis`
		* New permissible value: `Methicillin-Resistant Staphylococcus aureus (MRSA)`
		* New permissible value: `Staph Osteomyelitis`
		* New permissible value: `Thyroid Disease, Non-Cancer`
		* New permissible value: `Urinary Tract Infection`
	* Changes made to `diabetes_treatment_type`
		* New permissible value: `Linagliptin`
	* Changes made to `risk_factor`
		* New permissible value: `Adenomyosis`
		* New permissible value: `Human Herpesvirus-6 (HHV-6)`
		* New permissible value: `Human Herpesvirus-8 (HHV-8)`
	* Changes made to `undescended_testis_corrected_age`
		* New property: `maximum`
* Altered `molecular_test` Entity
	* New property: `hpv_strain`
	* Changes made to `antigen`
		* New permissible value: `TAG-72`
	* Changes made to `gene_symbol`
		* New permissible value: `TLR2`
	* Changes made to `laboratory_test`
		* New permissible value: `Erythrocyte Sedimentation Rate`
		* New permissible value: `Prothrombin Time`
	* Changes made to `second_gene_symbol`
		* New permissible value: `TLR2`
	* Changes made to `test_result`
		* New permissible value: `Stable`
	* Changes made to `test_units`
		* New permissible value: `cells/mL`
		* New permissible value: `mcg/L`
		* New permissible value: `mg/24 hr`
* Altered `somatic_annotation_workflow` Entity
	* Changes made to `workflow_type`
		* New permissible value: `GATK4 MuTect2 Tumor-Only Annotation`
* Altered `treatment` Entity
	* New property: `clinical_trial_indicator`
	* New property: `course_number`
	* New property: `drug_category`
	* New property: `embolic_agent`
	* New property: `lesions_treated_number`
	* New property: `number_of_fractions`
	* New property: `prescribed_dose`
	* New property: `protocol_identifier`
	* New property: `treatment_dose_max`
	* `route_of_administration`changed to array type.
	* Changes made to `treatment_dose_units`
		* New permissible value: `AUC`
		* New permissible value: `g/day`
		* New permissible value: `g/m2`
		* New permissible value: `IU/kg`
		* New permissible value: `IU/mg`
		* New permissible value: `mCi`
		* New permissible value: `mEq`
		* New permissible value: `mg/day`
		* New permissible value: `mg/dL`
		* New permissible value: `mg/kg`
		* New permissible value: `mg/kg/day`
		* New permissible value: `mg/m2`
		* New permissible value: `mg/m2/day`
		* New permissible value: `mg/m2/wk`
		* New permissible value: `mg/mL`
		* New permissible value: `mg/wk`
		* New permissible value: `mIU`
		* New permissible value: `mL`
		* New permissible value: `ug`
		* New permissible value: `ug/m2`
		* New permissible value: `Wafer`
	* Changes made to `treatment_type`
		* New permissible value: `Peptide Receptor Radionuclide Therapy (PRRT)`
	* Changes made to `therapeutic_agents`
		* New permissible value: `Canakinumab`
		* New permissible value: `Rivaroxaban`
* Altered `somatic_mutation_calling_workflow` Entity
	* Changes made to `workflow_type`
		* New permissible value: `GATK4 MuTect2 Tumor-Only`
* Altered `submitted_aligned_reads` Entity
	* Changes made to `experimental_strategy`
		* New permissible value: `m6A MeRIP-Seq`
		* Removed permissible value: `m6A RNA Methylation`

### Known Issues and Workarounds

* The `mitotic_count` field in `diagnosis` is erroneously set as "deprecated" and does not appear in the dictionary viewer.  This field can be uploaded successfully without issue and will appear in the dictionary viewer at a later release. <!--SV-1939-->

## v2.4.1

* __GDC Product__: GDC Data Dictionary
* __Release Date__: August 23, 2021

### New Features and Changes

* Altered `pathology_detail` Entity
	* Changes made to `additional_pathology_findings`
		* New permissible value: `Adenomyosis`
		* New permissible value: `Atrophic endometrium`
		* New permissible value: `Atypical hyperplasia/Endometrial intraepithelial neoplasia (EIN)`
		* New permissible value: `Autoimmune atrophic chronic gastritis`
		* New permissible value: `Asbestos bodies`
		* New permissible value: `Benign endocervical polyp`
		* New permissible value: `Bilateral ovaries with endometriotic cyst and surface adhesions`
		* New permissible value: `Carcinoma in situ`
		* New permissible value: `Cirrhosis`
		* New permissible value: `Clostridioides difficile (c. diff)`
		* New permissible value: `Colonization; bacterial`
		* New permissible value: `Colonization; fungal`
		* New permissible value: `Cyst(s)`
		* New permissible value: `Diffuse and early nodular diabetic glomerulosclerosis`
		* New permissible value: `Dysplasia; high grade`
		* New permissible value: `Dysplasia; low grade`
		* New permissible value: `Endometrial polyp`
		* New permissible value: `Endometriosis`
		* New permissible value: `Endometroid carcinoma with local mucinous differentiation`
		* New permissible value: `Endosalpingiosis`
		* New permissible value: `Epithelial dysplasia`
		* New permissible value: `Epithelial hyperplasia`
		* New permissible value: `Gallbladder adenomyomatosis`
		* New permissible value: `Glomerular disease`
		* New permissible value: `Hyperkeratosis`
		* New permissible value: `Inflammation`
		* New permissible value: `Intestinal metaplasia`
		* New permissible value: `Keratinizing dysplasia; mild`
		* New permissible value: `Keratinizing dysplasia; moderate`
		* New permissible value: `Keratinizing dysplasia; severe (carcinoma in situ)`
		* New permissible value: `Leiomyoma`
		* New permissible value: `Leiomyomata w/ degenerative changes`
		* New permissible value: `Nonkeratinizing dysplasia; mild`
		* New permissible value: `Nonkeratinizing dysplasia; moderate`
		* New permissible value: `Nonkeratinizing dysplasia; severe (carcinoma in situ)`
		* New permissible value: `Other`
		* New permissible value: `PD-L1 CPS (223C LDT) - 20%`
		* New permissible value: `Platinum-resistant`
		* New permissible value: `Pleural plaque`
		* New permissible value: `Pulmonary interstitial fibrosis`
		* New permissible value: `Sialadenitis`
		* New permissible value: `Sinonasal papilloma`
		* New permissible value: `Squamous metaplasia`
		* New permissible value: `Squamous papilloma; solitary`
		* New permissible value: `Squamous papillomatosis`
		* New permissible value: `Tubular (papillary) adenoma(s)`
		* New permissible value: `Tumor-associated lymphoid proliferation`
		* New permissible value: `Tumor has rough spikey edges`
		* Removed permissible value: `Pleurodesis, Talc`
		* Removed permissible value: `Pleurodesis, NOS`
* Altered `diagnosis` Entity
	* `primary_disease` field changed to `diagnosis_is_primary_disease` for clarity.
* Altered `treatment` Entity
	* Changes made to `treatment_type`
		* New permissible value: `Pleurodesis, Talc`
		* New permissible value: `Pleurodesis, NOS`
* Altered `molecular_test` Entity
	* Changes made to `test_units`
		* New permissible value: `mm^2`
		* New permissible value: `ng/mL`
		* New permissible value: `percent`
		* New permissible value: `count x10^9/L`

### Known Issues and Workarounds

* The `mitotic_count` field in `diagnosis` is erroneously set as "deprecated" and does not appear in the dictionary viewer.  This field can be uploaded successfully without issue and will appear in the dictionary viewer at a later release. <!--SV-1939-->


## v2.4.0

* __GDC Product__: GDC Data Dictionary
* __Release Date__: June 21, 2021

### New Features and Changes

* Added mix-max limitations on property values in `demographic`, `portion`, `aliquot`, `family_history`, `slide`, `follow_up`, `read_group`, `sample`, `analyte`, `exposure`, `diagnosis`, `treatment`, `molecular_test`
* Altered `submitted_unaligned_reads`, `submitted_aligned_reads`, `annotated_somatic_mutation`, `simple_somatic_mutation`, `masked_somatic_mutation`, `submitted_genomic_profile`, `aligned_reads`, `aggregated_somatic_mutation`, `simple_germline_variation` Entities
	* Changes made to `experimental_strategy`
		* Removed permissible value: `Low Pass WGS`
* Altered `follow_up` Entity
	* New property: `eye_color`
	* New property: `history_of_tumor`
	* New property: `history_of_tumor_type`
	* New property: `undescended_testis_corrected`
	* New property: `undescended_testis_corrected_age`
	* New property: `undescended_testis_corrected_laterality`
	* New property: `undescended_testis_corrected_method`
	* New property: `undescended_testis_history`
	* New property: `undescended_testis_history_laterality`
	* Changes made to `comorbidity`
		* New permissible value: `Dermatomyosis`
		* New permissible value: `Herpes Zoster`
		* New permissible value: `Varicella Zoster Virus`
	* Changes made to `risk_factor`
		* New permissible value: `Dermatomyosis`
		* New permissible value: `Herpes Zoster`
		* New permissible value: `Varicella Zoster Virus`
* Altered `read_group` Entity
	* Changes made to `instrument_model`
		* New permissible value: `Illumina NovaSeq 6000`
	* Changes made to `single_cell_library`
		* New permissible value: `Chromium scATAC v1 Library`
	* Changes made to `target_capture_kit`
		* New permissible value: `Twist Human Comprehensive Exome`
* Altered `somatic_aggregation_workflow` Entity
	* Added link: `simple_somatic_mutations`
	* Changes made to `workflow_type`
		* New permissible value: `CaVEMan Variant Aggregation and Masking`
* Altered `sample` Entity
	* Changes made to `sample_type`
		* New permissible value: `Next Generation Cancer Model Expanded Under Non-conforming Conditions`
	* Changes made to `sample_type_id`
		* New permissible value: `87`
* Altered `germline_mutation_calling_workflow` Entity
	* New link: `submitted_genotyping_arrays`
	* Changes made to `workflow_type`
		* New permissible value: `Birdseed`
* Altered `slide_image` Entity
	* Changes made to `data_format`
		* New permissible value: `JPEG 2000`
* Altered `exposure` Entity
	* New property: `exposure_duration_years`
	* New property: `parent_with_radiation_exposure`
	* Changes made to `exposure_type`
		* New permissible value: `Radiation`
* Altered `simple_germline_variation` Entity
	* New property: `platform`
	* Changes made to `data_format`
		* New permissible value: `TSV`
* Altered `pathology_detail` Entity
	* New property: `consistent_pathology_review`
	* New property: `residual_tumor`
	* New property: `size_extraocular_nodule`
	* New property: `tumor_thickness`
* Altered `diagnosis` Entity
	* New property: `adrenal_hormone`
	* New property: `primary_disease`
	* Changes made to `ajcc_pathologic_stage`
		* New permissible value: `Stage IIIA1`
		* New permissible value: `Stage IIIA2`
	* Changes made to `metastasis_at_diagnosis_site`
		* New permissible value: `Bladder`
		* New permissible value: `Bronchus`
		* New permissible value: `Head, Face or Neck, NOS`
		* New permissible value: `Lymph Node, Regional`
		* New permissible value: `Lymph Node, Subcarinal`
	* Changes made to `method_of_diagnosis`
		* New permissible value: `Exoresection`
	* Changes made to `morphology`
		* New permissible value: `8246/6`
		* New permissible value: `8380/6`
		* New permissible value: `8461/6`
		* New permissible value: `8522/6`
	* Changes made to `sites_of_involvement`
		* New permissible value: `Mesothelium`
* Altered `treatment` Entity
	* New property: `route_of_administration`
	* Changes made to `treatment_arm`
		* New permissible value: `A081801`
	* Changes made to `therapeutic_agents`
		* New permissible value: `Interferon Alfa-2B`
		* New permissible value: `Levoleucovorin Calcium`
		* New permissible value: `Mistletoe Extract`
* Altered `somatic_mutation_calling_workflow` Entity
	* Changes made to `workflow_type`
		* New permissible value: `Strelka2 RNA`
* Altered `molecular_test` Entity
	* New property: `days_to_test`
	* New link: `diagnoses`
	* Changes made to `antigen`
		* New permissible value: `FMC-7`
		* New permissible value: `Kappa, Surface`
		* New permissible value: `Lambda, Surface`
	* Changes made to `gene_symbol`
		* New permissible value: `AQP1`
		* New permissible value: `CALB2`
		* New permissible value: `DNTT`
		* New permissible value: `EPCAM`
		* New permissible value: `GCET1`
		* New permissible value: `PDPN`
		* New permissible value: `PTGS2`
	* Changes made to `laboratory_test`
		* New permissible value: `BG8`
		* New permissible value: `Circulating Endothelial Cells`
		* New permissible value: `Cytokeratin 5`
		* New permissible value: `Cytokeratin 6`
		* New permissible value: `Dopamine-Secreting`
		* New permissible value: `Epinephrine-Secreting`
		* New permissible value: `Metanephrine-Secreting`
		* New permissible value: `Methoxytyramine-Secreting`
		* New permissible value: `Microsatellite Instability`
		* New permissible value: `Norepinephrine-Secreting`
		* New permissible value: `Normetanephrine-Secreting`
		* New permissible value: `Serum Mesothelin`
		* New permissible value: `TAG-72`

### Known Issues and Workarounds

* The `mitotic_count` field in `diagnosis` is erroneously set as "deprecated" and does not appear in the dictionary viewer.  This field can be uploaded successfully without issue and will appear in the dictionary viewer at a later release. <!--SV-1939-->


## v2.3.0

* __GDC Product__: GDC Data Dictionary
* __Release Date__: January 5, 2021

### New Features and Changes
* Altered `submitted_unaligned_reads` Entity
	* Changes made to `experimental_strategy`
		* New permissible value: `HiChIP`
		* New permissible value: `m6A RNA Methylation`
		* New permissible value: `scATAC-Seq`
	* Changes made to `read_pair_number`
		* New permissible value: `R3`
* Altered `submitted_aligned_reads` Entity
	* Changes made to `experimental_strategy`
		* New permissible value: `HiChIP`
		* New permissible value: `m6A RNA Methylation`
		* New permissible value: `scATAC-Seq`
* Altered `follow_up` Entity
	* Changes made to `comorbidity`
		* New permissible value: `Abnormal Glucose Level`
		* New permissible value: `Chronic Fatigue Syndrome`
		* New permissible value: `Clonal Hematopoiesis`
		* New permissible value: `Fibromyalgia`
		* New permissible value: `Gastritis`
	* Changes made to `risk_factor`
		* New permissible value: `Abnormal Glucose Level`
		* New permissible value: `Chronic Kidney Disease`
		* New permissible value: `Escherichia coli`
		* New permissible value: `Gastritis`
		* New permissible value: `Skin Rash`
* New Entity: `masked_methylation_array`
* Altered `read_group` Entity
	* New property: `chipseq_antibody`
	* New property: `fragmentation_enzyme`
	* Removed property: `RIN`
	* Changes made to `library_strategy`
		* New permissible value: `HiChIP`
		* New permissible value: `m6A RNA Methylation`
		* New permissible value: `scATAC-Seq`
* Altered `sample` Entity
	* New property: `sample_ordinal`
	* Changes made to `composition`
		* New permissible value: `Mixed Adherent Suspension`
* Altered `analyte` Entity
	* New property: `experimental_protocol_type`
	* New property: `rna_integrity_number`
* Altered `pathology_detail` Entity
	* New property: `additional_pathology_findings`
	* New property: `necrosis_percent`
	* New property: `necrosis_present`
	* New property: `rhabdoid_percent`
	* New property: `rhabdoid_present`
	* New property: `sarcomatoid_percent`
	* New property: `sarcomatoid_present`
	* Changes made to `dysplasia_degree`
		* New permissible value: `Mild`
		* New permissible value: `Moderate`
		* New permissible value: `Severe`
	* Changes made to `dysplasia_type`
		* New permissible value: `Epithelial`
		* New permissible value: `Keratinizing`
		* New permissible value: `Nonkeratinizing`
	* Changes made to `lymph_node_involvement`
* Altered `diagnosis` Entity
	* New property: `ann_arbor_b_symptoms_described`
	* Changes made to `ajcc_clinical_stage`
		* New permissible value: `Stage IA3`
	* Changes made to `ajcc_pathologic_m`
		* New permissible value: `M1d`
	* Changes made to `metastasis_at_diagnosis_site`
		* New permissible value: `Gastrointestinal Tract`
		* New permissible value: `Heart`
		* New permissible value: `Neck`
		* New permissible value: `Retroperitoneum`
		* New permissible value: `Urethra`
		* New permissible value: `Uterine Adnexa`
		* New permissible value: `Vertebral Canal`
		* New permissible value: `Vulva, NOS`
	* Changes made to `morphology`
		* New permissible value: `8249/6`
		* New permissible value: `8800/6`
	* Changes made to `supratentorial_localization`
		* New permissible value: `Frontal lobe`
		* New permissible value: `Occipital lobe`
		* New permissible value: `Parietal lobe`
		* New permissible value: `Temporal lobe`
* Altered `treatment` Entity
	* Changes made to `treatment_dose_units`
		* New permissible value: `mg`
	* Changes made to `treatment_type`
		* New permissible value: `Radiation, Hypofractionated`
		* New permissible value: `Radiation, Mixed Photon Beam`
		* New permissible value: `Radiation, Photon Beam`
	* Changes made to `therapeutic_agents`
		* Updated enum list.
* Altered `case` Entity
	* Updated `disease_type` and `primary_site` enum values.
* Altered `rna_expression_workflow` Entity
	* Changes made to `workflow_type`
		* New permissible value: `STAR - Smart-Seq2 Raw Counts`
		* New permissible value: `STAR - Smart-Seq2 Filtered Counts`
* Altered `molecular_test` Entity
	* New link: `slides`
	* Changes made to `antigen`
		* New permissible value: `Ki67`
	* Changes made to `gene_symbol`
		* New permissible value: `CHGA`
		* New permissible value: `SYP`
	* Changes made to `molecular_consequence`
		* New permissible value: `Exon Variant`
	* Changes made to `second_gene_symbol`
		* New permissible value: `CHGA`
		* New permissible value: `SYP`
* Altered `aligned_reads` Entity
	* Changes made to `experimental_strategy`
		* New permissible value: `HiChIP`
		* New permissible value: `scATAC-Seq`

### Bugs Fixed Since Last Release
* None

## v2.2.0

* __GDC Product__: GDC Data Dictionary
* __Release Date__: July 2, 2020

### New Features and Changes
* Added NCIt codes associated with enumeration values in `diagnosis` entity type
* Added `pathology_detail` entity <!--DICT-16-->
* Modified `treatment` entity <!--DICT-1-->
    - Added new enumerations for `therapeutic_agents` property
        * `Itraconazole`
        * `Tipiracil`
        * `Tipiracil Hydrochloride`
        * `Zirconium Zr 89 Panitumumab`
        * `Estradiol mustard`
        * `Progestational IUD`
        * `PD-1 Inhibitor`
        * `IGF-1R Inhibitor`
        * `CDK4/6 Inhibitor`
        * `ALK Inhibitor`
* Modified `annotation` entity <!--DICT-2 -->
    - Added `many_to_many` link to `copy_number_estimate` entity
* Modified `diagnosis` entity <!--DICT-6-->
    - Added new properties:
        * `eln_risk_classification`
        * `satellite_nodule_present`
        * `who_cns_grade`
        * `who_nte_grade`
        *  `sites_of_involvement`
    - Added new permissible values to `ajcc_pathologic_stage` property
        * `Stage IA3`
    - Added new permissible values to `classification_of_tumor` property
        * `Progression`
    - Added new permissible values to `metastasis_at_diagnosis_site` property
        * `Esophagus`
    - Added new permissible values to `morphology` property
    - Removed properties as `required`
        * `days_to_last_follow_up`
        * `tumor_grade`
        * `progression_or_recurrence`
        * `days_to_recurrence`
        * `days_to_last_known_disease_status`
        * `last_known_disease_status`
    - Deprecated properties
        * `mitotic_count`
        * `papillary_renal_cell_type`
        * `micropapillary_features`
        * `non_nodal_regional_disease`
        * `non_nodal_tumor_deposits`
* Modified `sample` entity <!--DICT-7-->
    - Changed description of properties
        * `days_to_collection`
        * `days_to_sample_procurement`
    - Deprecated properties
        * `is_ffpe`
        * `oct_embedded`
    - Added new permissible values to `sample_type` property
        * `Mixed Adherent Suspension Saliva`
* Modified `follow_up` entity <!--DICT-11-->
    - Added new properties
        * `procedures_performed`
        * `hormonal_contraceptive_type`
        * `hormone_replacement_therapy_type`
    - Added new permissible values to `comorbidity` property
    - Added new permissible values to `risk_factor` property
    - Added new permissible values to `evidence_of_recurrence_type` property
    - Added new permissible values to `aids_risk_factors` property
* Modified `exposure` entity <!--DICT-12-->
    - Added new properties
        * `smokeless_tobacco_quit_age`
        * `alcohol_type`
    - Added new permissible values to `exposure_type` property
        * `Wood Dust`
        * `Smoke`
    - Added new permissible value to `type_of_smoke_exposure` property
        *  `Tobacco smoke, NOS`
* Modified `submitted_unaligned_reads` property <!--DICT-13-->
    - Removed permissible value from `read_pair_number` property
        * `I1`
* Modified `aliquot` entity <!--DICT-14-->
    - Added new permissible value to `analyte_type` property
        * `Nuclei RNA`
* Modified `rna_expression_workflow` entity <!--DICT-15-->
    - Removed permissible value from `workflow_type` property
        * `STAR - Smart-Seq2 Counts`
    - Added new permissible values to `workflow_type` property
        * `STAR - Smart-Seq2 Gene Counts`
        * `STAR - Smart-Seq2 GeneFull Counts`
* Modified `molecular_test` entity <!--DICT-17-->
    - Added new properties
        * `mitotic_count`
        * `mitotic_total_area`
        * `biospecimen_volume`
    - Added new permissible values to `gene_symbol` property
    - Added new permissible values to `second_gene_symbol` property
    - Added new permissible values to `antigen` property
    - Added new permissible values to `laboratory_test` property 	  

### Bugs Fixed Since Last Release
* None

## v2.1.0

* __GDC Product__: GDC Data Dictionary
* __Release Date__: March 10, 2020

### New Features and Changes

* Added NCIt codes associated with enumeration values in `exposure`, `family_history`, and `demographic` entity types
* Restructured dictionary for consistency across types
* __read_group entity__
* Added 5 new target capture kits <!--DAT-2755-->
    - `Custom Twist Broad PanCancer Panel - 396 Genes`
    - `Nextera DNA Exome`
    - `Custom Twist Broad Exome v1.0 - 35.0 Mb`
    - `Custom SureSelect CGCI-HTMCP-CC KMT2D And Hotspot Panel - 37.0 Kb`
    - `TruSeq RNA Exome`
* Add one new enum to `library_strategy` <!--DAT-2793-->
    - `scRNA-Seq`
* __structural_variation entity__ <!--DAT-2703-->
* Added `VCF` to `data_format` property
* Added links to `structural_variation` from `somatic_mutation_index` <!--DAT-2704-->
* Added links to `aligned_reads` from `somatic_copy_number` workflow <!--DAT-2705-->
* __copy_number_segment entity__
* Added new permissible values to `experimental_strategy` property <!--DAT-2706-->
    - `WGS`
    - `WXS`
* Added new enum to `experimental_strategy` <!--DAT-2788-->
    - `WGS`
* __aligned_reads entity__
* Added 2 new properties <!--DAT-2787-->
    - `tumor_ploidy`
    - `tumor_purity`
* __treatment entity__
* Added enumeration to `therapeutic_agent` property <!--DAT-2727-->
* __copy_number_estimate entity__ <!--DAT-2729-->
* Added new enum to `data_format`
    - `TSV`
* __demographic entity__	 <!--DAT-2733-->
* Added new property `country_of_residence_at_enrollment`
* __family_history entity__ <!--DAT-2738-->
* Added new permissible values to `relationship_primary_diagnosis` property
*  __follow_up entity__	 <!--DAT-2739-->
* Added new properties
    - `body_surface_area`
    - `recist_targeted_regions_number`
    - `recist_targeted_regions_sum`
    - `adverse_event_grade`
    - `cd4_count`
    - `imaging_type`
    - `scan_tracer_used`
    - `nadir_cd4_count`
    - `hiv_viral_load`
    - `aids_risk_factors`
    - `haart_treatment_indicator`
    - `immunosuppressive_treatment_type`
    - `evidence_of_recurrence_type`
    - `imaging_result`
    - `hormonal_contraceptive_use`
    - `pregnancy_outcome`
    - `hysterectomy_type`
    - `hysterectomy_margins_involved`
    - `days_to_imaging`
    - `cdc_hiv_risk_factors`
    - `risk_factor`
* Added new enum to `days_to_follow_up`
    - `null`
* __molecular_test entity__ <!--DAT-2741-->
* Added new permissible values to various properties
    - `laboratory_test`
    - `second_gene_symbol`
    - `molecular_consequence`
    - `biospecimen_type`
    - `molecular_analysis_method`
    - `gene_symbol`
    - `clonality`
* __sample entity__ <!--DAT-2742-->
* Added new permissible values to various properties
    - `method_of_sample_procurement`
    - `biospecimen_anatomic_site`
    - `tumor_descriptor`
* Added new property
    - `tissue_collection_type`
* __treatment entity__ <!--DAT-2745-->
* Added new properties
    - `treatment_arm`
    - `reason_treatment_ended`
    - `number_of_cycles`
    - `treatment_effect_indicator`
    - `treatment_dose`
    - `treatment_dose_units`
    - `treatment_frequency`
    - `chemo_concurent_to_radiation`
* Added new permissible values to various properties
    - `therapeutic_agent`
    - `treatment_effect`
    - `treatment_intent_type`
* __slide entity__ <!--DAT-2747-->
* Added new properties
    - `percent_sarcomatoid_features`
    - `percent_rhabdoid_features`
    - `prostatic_chips_total_count`
    - `prostatic_chips_positive_count`
    - `prostatic_involvement_percent`
    - `bone_marrow_malignant_cells`
    - `percent_follicular_component`
    - `tissue_microarray_coordinates`
* __diagnosis entity__ <!--DAT-2735-->
* Added new properties
    - `tumor_depth`
    - `margin_distance`
    - `transglottic_extension`
    - `margins_involved_site`
    - `gleason_grade_tertiary`
    - `papillary_renal_cell_type`
    - `gleason_patterns_percent`
    - `greatest_tumor_dimension`
    - `lymph_node_involved_site`
    - `pregnant_at_diagnosis`
    - `figo_staging_edition_year`
* Added new permissible values to various properties
    - `classification_of_tumor`
    - `figo_stage`
    - `tumor_grade`
    - `metastasis_at_diagnosis_site`
    - `gleason_grade_group`
    - `tissue_or_organ_of_origin`
    - `morphology`
* __exposure entity__ <!--DAT-2737-->
* Added new properties
    - `secondhand_smoke_as_child`
    - `exposure_type`
    - `type_of_tobacco_used`
    - `exposure_duration`
    - `tobacco_use_per_day`
    - `age_at_onset`
    - `marijuana_use_per_week`
    - `tobacco_smoking_status`
* __case entity__ <!--DAT-2746-->
* Added new properties
    - `consent_type`
    - `days_to_consent`
* Added new permissible values to `index_date`
* Added `scRNA-Seq` as new enum to `experimental_strategy` in 4 entities <!--DAT-2794-->
    - `submitted_unaligned_reads`
    - `submitted_aligned_reads`
    - `aligned_reads`
    - `gene_expression`
* Added 3 new permissible values to `workflow_type` in alignment_workflow <!--2795-->
    - `CellRanger - 10x Chromium`
    - `STAR - Smart-Seq2`
    - `zUMIs - Smart-Seq2`
* Added 4 new permissible values to `workflow_type` in rna_expression_workflow <!--DAT-2796-->
    - `CellRanger - 10x Raw Counts`
    - `CellRanger - 10x Filtered Counts`
    - `STAR - Smart-Seq2 Counts`
    - `zUMIs - Smart-Seq2 Counts`
* Add one new integer property to `read_group` <!--DAT-2798-->
    -  `number_expect_cells`
* Add one new enum to `read_pair_number` in `submitted_unaligned_reads` <!--DAT-2799-->
    - `I1`
* Add one new enum property in `read_group` node <!--DAT-2800-->
    - `Chromium 3' Gene Expression v2 Library`
    - `Chromium 3' Gene Expression v3 Library`
    - `Smart-Seq2`
* Created new node `expression_analysis_workflow` <!--DAT-2802-->
* Created new node `secondary_expression_analysis` <!--2803-->
* Add one new enum to `data_format` in `gene_expression` <!--DAT-2804-->
    -  `MEX`

### Bugs Fixed Since Last Release
* None

## v2.0.0 <!--REQ-393-->

* __GDC Product__: GDC Data Dictionary
* __Release Date__: January 30, 2020

### New Features and Changes
* The API that includes the GDC data dictionary now uses Python 3.

### Bugs Fixed Since Last Release

* None


## v.1.18.1

* __GDC Product__: GDC Data Dictionary
* __Release Date__: November 6, 2019

### New Features and Changes
* Added new permissible value `deleted` to property `file_state`

### Bugs Fixed Since Last Release

* None


## v.1.18

* __GDC Product__: GDC Data Dictionary
* __Release Date__: July 31, 2019

### New Features and Changes

* Added new entities <!--DAT-2450-->
    - `protein_expression_quantification`  <!--DAT-2237-->
    - `submitted_genotyping_array`  <!--DAT-2327-->
    - `somatic_copy_number_workflow` <!--DAT-2330-->
* Add links in data model to `somatic_copy_number_workflow` from `copy_number_segment`  <!--DAT-2331-->
* Add links in data model to `somatic_copy_number_workflow` from `copy_number_estimate`	 <!--DAT-2415-->
* Add links in data model from `annotated_somatic_mutation` to `genomic_profile_harmonization_workflow` <!--DAT-2333-->
* Modified `copy_number_segment` entity
    - Add new data type <!--DAT-2332-->
        - `Allele-specific Copy Number Segment`
* Update data dictionary to support new annotation classifications <!--DAT-2238-->
* Fixed typo in `sample` entity schema <!--DAT-2247-->
* Unrequired project link for aliquot level MAFs  <!--DAT-2260-->
* Added NCIt codes for gender values  <!--DAT-2314-->
* Changed description of `masked_somatic_mutation` and `aggregated_somatic_mutation` nodes to be the same  <!--DAT-2311-->
* Modified `archive` entity <!--DAT-2259-->
    - Set `downloadable` property to `true`
* Modified `publication`  entity<!--DAT-2259-->
    - Set `downloadable` property to `true`
* Modified `filtered_copy_number_segment`  entity<!--DAT-2259-->
    - Set `downloadable` property to `true`
* Modified `aligned_reads` entity
    - Added `MSI properties` as new property
* Modified `read_group` entity
    - Added `Custom SureSelect Human All Exon v1.1 Plus 3 Boosters` as new permissible value for `target_capture_kit` field <!--DAT-2003-->
    - Added `SeqCap EZ Human Exome v3.0` as new permissible value for `target_capture_kit` field <!--DAT-2310-->
    - Added new permissible values for `instrument_model`
        - `Unknown` <!--DAT-2411-->
        - `Not Reported` <!--DAT-2411-->
        - `Ion Torrent S5` <!--DAT-2338-->
* Modified `biospecimen_supplement` entity <!--  DAT-2266-->
    - Added `CDC JSON` as new permissible data format
* Modified `demographic` entity <!--DAT-2303-->
    - Added new property
        - `age_is_obfuscated`
* Modified `demographic` entity <!--DAT-2303-->
    - Added new properties
        - `cause_of_death_source` <!--DAT-2341-->
        - `occupation_duration_years` <!--DAT-2342-->
* Modified `diagnosis` entity
    - Added new properties
        - `non_nodal_regional_disease` <!--DAT-2353-->
        - `non_nodal_tumor_deposits` <!--DAT-2354-->
        - `ovarian_specimen_status` <!--DAT-2356-->
        - `ovarian_surface_involvement` <!--DAT-2357-->
        - `percent_tumor_invasion` <!--DAT-2360-->
        - `peritoneal_fluid_cytological_status` <!-- DAT-2361-->      
        - `breslow_thickness` <!--DAT-2345-->
        - `international_prognostic_index` <!--DAT-2349-->
        - `largest_extrapelvic_peritoneal_focus` <!--DAT-2350-->
        - `mitotic_count` <!--DAT-2352-->     
    - Removed permissible values from `primary_diagnosis`  <!--DAT-2377 -->
    - Removed permissible values from `site_of_resection_or_biopsy` <!--DAT-2378-->
    - Removed `tumor_stage` as required property <!--DAT-2397-->
    - Added new permissible values for
        - `ajcc_pathologic_stage` <!--DAT-2368-->
        - `metastasis_at_diagnosis_site` <!--DAT-2370-->
    - Migrated values not part of permissible values to `Not Reported` for following properties
          - `primary_diagnosis` <!--DAT-2412-->
          - `site_of_resection_or_biopsy` <!--DAT-2413 -->
          - `tumor_grade` <!--DAT-2414-->
          - `tissue_or_organ_of_origin` <!--DAT-2438-->
    - Removed permissible values from
          - `tissue_or_organ_of_origin` <!--DAT-2426-->
          - `morphology` <!--  DAT-2427-->
- Modified `structural_variant_calling_workflow` entity <!--DAT-233-->
    - Added new `workflow_type`
- Modified `structural_variant` entity
    - Added `BEDPE` to `data_format` as permissible value <!--DAT-2336-->
- Modified `molecular_test` entity
    - Added new permissible values for `gene_symbol` <!--DAT-2373-->
    - Added new permissible values for `test_result`<!--DAT-2374-->
    - Added new permissible values for `antigen`<!--DAT-2372-->
    - Added new property `pathogenicity`  <!--DAT-2366-->
- Modified `follow_up` entity
    - Added new permissible value for `risk_factor` <!--DAT-2371-->
- Modified `sample` entity
    - Added new permissible values for `method_of_sample_procurement` <!--DAT-2376-->
- Modified `genomic_profile_harmonization_workflow` entity
    - Added new `workflow_type` permissible values <!--DAT-2388-->
- Modified `somatic_mutation_calling_workflow` entity
    - Added new `workflow_type` permissible values <!--DAT-2389-->
    - Modified `somatic_annotation_workflow` entity
    - Added new `workflow_type` permissible values <!--DAT-2391-->
- Modified `case` entity
    - Added new permissible value to `disease_type` <!--DAT-2390-->
        - `Not Applicable`
- Modified `copy_number_estimate` entity
    - Added new permissible value to `experimental_strategy` <!--DAT-2334-->
        -  `WXS`
- Modified `family_history` entity
    - Added new property <!--DAT-2365-->
        - `relatives_with_cancer_history_count`
- Modified `sample` entity <!--DAT-2418-->
    - Added new permissible values
        - `sample_type`
        - `sample_type_id`

### Bugs Fixed Since Last Release

* None

## v.1.17

* __GDC Product__: GDC Data Dictionary
* __Release Date__: June 5, 2019


### New Features and Changes

* Deleted vital status, days_to_birth, and days_to_death from Diagnosis node.  Data submission and data requests should all be directed to the corresponding properties on the Demographic Node.

### Bugs Fixed Since Last Release

* None



## v.1.16

* __GDC Product__: GDC Data Dictionary
* __Release Date__: April 17, 2019


### New Features and Changes

* Updates to the Data Dictionary Search Tool
* Added new bioinformatics workflow for methylation arrays (Sesame) <!--DAT-2025-->
* Changed `somatic_mutation_calling_workflow` link from `one_to_many` to `many_to_many`	<!--DAT-1697-->
* Modified `read_group` entity
    - Added `SeqCap EZ Human Exome v2.0` as new permissible value for `target_capture_kit` field <!--DAT-1827-->
    - Added `Custom SureSelect Human All Exon v1.1 Plus 3 Boosters` as new permissible value `target_capture_kit` field <!--DAT-2002-->
    - Added `Custom SureSelect CGCI-HTMCP-CC Panel - 19.7 Mb`  as new permissible value `target_capture_kit` field <!--DAT-2013-->
* Modified `case` entity
    - Updated the description for the `primary_site` field <!--DAT-1928-->
    - Added new permissible value to `lost_to_followup` field <!--DAT-2133-->
* Modified `molecular_test` entity
    - Removed properties with genomic coordinates <!--DAT-1991-->
    - Add new permissible values to `test_result`	<!--TT-918-->
    - Added `second_exon` as new property<!--TT-919-->
* Modified `aligned_reads_index` entity
    - Made these files not submittable <!--DAT-1985-->
* Modified `somatic_mutation_index` entity
    - Made these files not submittable <!--DAT-1986-->
* Modified `sample` entity
    - Added new permissible values for `sample_type` <!--DAT-2017-->
        - `Blood Derived Cancer - Bone Marrow`
        - `Blood Derived Cancer - Peripheral Blood`
    - Added new permissible values to `sample_type_id` <!--DAT-2116-->
* Modified `diagnosis` entity
    - Added 6 new staging and grading properties for TCGA <!--DAT-2056-->
        - `igcccg_stage`
        - `masaoka_stage`
        - `gleason_grade_group`
        - `primary_gleason_grade`
        - `secondary_gleason_grade`
        - `weiss_assessment_score`
    - Made `vital_status` an optional field <!--DAT-2157-->
    - Removed deprecated properties <!--TT-939--><!--TT-938-->
        - `days_to_death`
        - `days_to_birth`
        - `cause_of_death`
        - `hiv_positive`
        - `days_to_hiv_diagnosis`
        - `ldh_normal_range_upper`
        - `new_event_type`
        - `hpv_status`
        - `hpv_positive_type`
        - `colon_polyps_history`
        - `progression_free_survival`
        - `progression_free_survival_event`
        - `overall_survival`
        - `days_to_treatment`
        - `ldh_level_at_diagnosis`
    - Added `vital_status` property to deprecated list<!--DAT-2158-->
* Modified `somatic_aggregation_workflow` entity
    - Added `Aliquot Ensemble Somatic Variant Merging and Masking` as new permissible value to `workflow_type`
* Modified `slide` entity
    - Updated the description for the `magnification` field <!--DAT-2125-->
* Modified `aliquot` entity
    - Updated the the description for several fields <!--DAT-2126-->
        - `selected_normal_low_pass_wgs`
        - `selected_normal_targeted_sequencing`
        - `selected_normal_wgs`
        - `selected_normal_wxs`
* Modified `follow-up` entity
    - Added new field `days_to_progression_free` <!--DAT-2135-->
* Modified `demographic` entity
    - Made `vital_status` a required field <!--DAT-2157-->
* Modified `exposure` entity <!--TT-926-->
    - Added new properties
        - `environmental_tobacco_smoke_exposure`
        - `respirable_crystalline_silica_exposure`
        - `coal_dust_exposure`
        - `type_of_smoke_exposure`
        - `type_of_tobacco_used`
        - `smoking_frequency`
        - `time_between_waking_and_first_smoke`
    - Removed `cigarettes_per_day` property from deprecated list<!--DAT-2062-->
* Modified `annotation` entity
    - Modified permissible values to `status` <!--TT-930-->
        - Approved
        - Rescinded

### Bugs Fixed Since Last Release

* None


## v.1.15

* __GDC Product__: GDC Data Dictionary
* __Release Date__: December 18, 2018


### New Features and Changes

* Removed `Raw Sequencing Data` and `Sequencing Data` as permissible values from `submitted_aligned_reads`, `submitted_unaligned_reads`, and `aligned_reads` <!--DAT-42--> <!--DAT-1904-->
* Deleted `aligned_reads_metrics` entity <!--DAT-1754-->
* Created new `raw_methylation_array` entity <!--DAT-1854-->
* Add regex validation to property `md5sum` for following entities: <!--DAT-1899-->
    - `slide_image`
    - `analysis_metadata`
    - `clinical_supplement`
    - `experiment_metadata`
    - `pathology_report`
    - `run_metadata`
    - `biospecimen_supplement`
    - `submitted_aligned_reads`
    - `submitted_genomic_profile`
    - `submitted_methylation_beta_value`
    - `submitted_tangent_copy_number`
    - `submitted_unaligned_reads`
* Modified `molecular_test` entity
    - Migrated data from `blood_test` to `laboratory_test` and `biospecimen_type` for all entities<!--TT-754-->
    - Added new property `intron` <!--DAT-1847-->
    - Deleted `blood_test` entity <!--DAT-1639-->
    - Added new permissible values for `gene_symbol`<!--DAT-1553-->
    - Added new permissible values for `antigen`<!--DAT-1662-->
    - Added new permissible values for `molecular_analysis_method` <!--DAT-1663-->
    - Added new permissible values for `variant_type` <!--DAT-1664-->
    - Added new permissible values for `test_result` <!--DAT-1665-->
    - Added new permissible values for `molecular_consequence` <!--DAT-1666-->
    - Added regex validation to property `transcript` <!--DAT-1916-->
    - Added regex validation to property `locus` <!--DAT-1874-->
    - Changed data type of `exon` property to be `string` with regex validation <!--DAT-1890-->
* Modified `diagnosis` entity
    - Added new fields
      - `tumor_focality`<!-- DAT-1832-->
      - `tumor_regression_grade` <!--DAT-1833-->
      - `lymph_nodes_tested` <!--DAT-1834-->
    - Added new permissible value for `primary_diagnosis` field<!--DAT-1879-->
    - Added min and max values to time-based properties <!--DAT-1885-->
    - Added new permissible value for `morphology` field <!--TT-818-->
* Modified `follow_up` entity
    - Added new permissible values for `ecog_performance_status`<!--DAT-1684-->
    - Added new permissible values for `comorbidity` <!--DAT-1766-->
    - Added new permissible values for `disease_response`<!--DAT-1840-->
    - Added new permissible values for `risk_factor`<!-- DAT-1841-->
    - Added min and max values to time-based properties <!--DAT-1884-->
    - Added new property:
      - `hepatitis_sustained_virological_response` <!--DAT-1845-->
    - Updated CDE, CDE version, description and URL for `comorbidity`<!--DAT-1911-->
    - Added a CDE for `days_to_comorbidity` <!--DAT-1912-->
    - Removed `reflux_treatment` property <!--DAT-1913-->
    - Add a new property:<!--DAT-1843-->
      - `risk_factor_treatment`
* Modified `aligned_reads` entity
  - Added new contamination properties <!--DAT-1749-->
    - `contamination`
    - `contamination_error`
* Modified `read_group` entity
  - Added new permissible values for `target_capture_kit` <!--DAT-1757-->
  - Updated description for property `instrument_model` <!--DAT-1763-->
  - Added new permissible values for `target_capture_kit` <!--DAT-1799-->
  - Added new permissible values for `library_strategy`<!--DAT-1814-->
  - Added regex validation to property `adapter_sequence` <!--DAT-1895-->
  - Added regex validation to property `multiplex_barcode`<!--DAT-1897-->
  - Allow users to enter null for property `read_length` <!--DAT-1908-->
  - Allow users to enter null for property `is_paired_end` <!--DAT-1909-->
* Modified `family_history` entity
  - Added new permissible values for `relationship_primary_diagnosis` <!--DAT-1765-->
  - Added min and max values to properties <!--DAT-1887-->
* Modified `case` entity
  - Add min and max values to properties <!--DAT-1888-->
  - Delete permissible value from `primary_site` <!-- DAT-1772-->
    - `Unknown Primary Site`
* Modified `analyte` entity
  - Corrected the description for fields `analyte_volume` to include microliters as unit <!--DAT-1801-->
* Modified `exposure` entity
  - Added new properties
    - `asbestos_exposure` <!-- DAT-1836-->
    - `radon_exposure` <!-- DAT-1837-->
* Modified `sample` entity
  - Added new permissible values to `method_of_sample_procurement` <!--DAT-1849-->
  - Added regex validation to `pathology_report_uuid` <!--DAT-1893-->
  - Change type from string to number for properties:
    - `intermediate_dimension` <!--DAT-1861-->
    - `longest_dimension` <!--DAT-1863-->
    - `shortest_dimension` <!--DAT-1865-->
    - `time_between_clamping_and_freezing` <!--DAT-1867-->
    - `time_between_excision_and_freezing` <!--DAT-1869-->
  - Add min and max to properties on the sample node <!--DAT-1883-->
  - Populated sample nodes that have no value for `tissue_type` to "Not Reported" <!--TT-658-->
* Modified `treatment` entity
  - Added a new property
    - `prior_treatment_effect` <!--DAT-1850-->
  - Add min and max values to properties <!--DAT-1889-->
* Modified `aliquot` entity
  - Corrected the description for fields `analyte_volume` to include microliters as unit <!--DAT-1859-->
* Modified `demographic` entity
  - Added min and max to properties <!--DAT-1886-->


### Bugs Fixed Since Last Release

* Fixed value of `pathology_report_uuid` on sample entity `7b29b034-86e4-4266-8657-036e96e04430` to satisfy regex requirements <!--DAT-1939-->
* Migrated a few unsupported values for sample.pathology_report_uuid, read_group.adapter_sequence, read_group.multiplex_barcode <!--DAT-1941-->

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
    - Added new permissible values for `antigen` field <!--DAT-1662-->
    - Added new permissible values to `molecular_analysis_method` <!--DAT-1663-->
    - Added new permissible values for `variant_type` field <!--DAT-1664-->
    - Added new permissible values to `test_result` <!--DAT-1665-->
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

* None

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
