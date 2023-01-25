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