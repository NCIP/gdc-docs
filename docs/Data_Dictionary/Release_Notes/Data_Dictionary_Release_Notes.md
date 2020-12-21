# Data Dictionary Release Notes

| Version | Date |
|---|---|
| [v.2.3.0](Data_Dictionary_Release_Notes.md#v220) | December 21, 2020 |
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

## v2.3.0

* __GDC Product__: GDC Data Dictionary
* __Release Date__: December 21, 2020

### New Features and Changes
* Altered `submitted_unaligned_reads` Entity
	* Changes made to `properties`
		* Changes made to `experimental_strategy`
			* New permissable value: `HiChIP`
			* New permissable value: `m6A RNA Methylation`
			* New permissable value: `scATAC-Seq`
		* Changes made to `read_pair_number`
			* New permissable value: `R3`
* Altered `annotation` Entity
	* Changes made to `links`
		* `masked_methylation_arrays` added to subgroup
		* `molecular_tests` added to subgroup
		* `pathology_details` added to subgroup
	* Changes made to `properties`
		* New property: `masked_methylation_arrays`
		* New property: `molecular_tests`
		* New property: `pathology_details`
* Altered `submitted_aligned_reads` Entity
	* Changes made to `properties`
		* Changes made to `experimental_strategy`
			* New permissable value: `HiChIP`
			* New permissable value: `m6A RNA Methylation`
			* New permissable value: `scATAC-Seq`
* Altered `follow_up` Entity
	* Changes made to `properties`
		* Changes made to `comorbidity`
			* New permissable value: `Abnormal Glucose Level`
			* New permissable value: `Chronic Fatigue Syndrome`
			* New permissable value: `Clonal Hematopoiesis`
			* New permissable value: `Fibromyalgia`
			* New permissable value: `Gastritis`
		* Changes made to `risk_factor`
			* New permissable value: `Abnormal Glucose Level`
			* New permissable value: `Chronic Kidney Disease`
			* New permissable value: `Escherichia coli`
			* New permissable value: `Gastritis`
			* New permissable value: `Skin Rash`
* New Entity: `masked_methylation_array`
* Altered `read_group` Entity
	* Changes made to `properties`
		* New property: `chipseq_antibody`
		* New property: `fragmentation_enzyme`
		* Removed property: `RIN`
		* Changes made to `library_strategy`
			* New permissable value: `HiChIP`
			* New permissable value: `m6A RNA Methylation`
			* New permissable value: `scATAC-Seq`
* Altered `sample` Entity
	* Changes made to `properties`
		* New property: `sample_ordinal`
		* Changes made to `composition`
			* New permissable value: `Mixed Adherent Suspension`
* Altered `analyte` Entity
	* Changes made to `properties`
		* New property: `experimental_protocol_type`
		* New property: `rna_integrity_number`
* Altered `pathology_detail` Entity
	* Changes made to `properties`
		* New property: `additional_pathology_findings`
		* New property: `necrosis_percent`
		* New property: `necrosis_present`
		* New property: `rhabdoid_percent`
		* New property: `rhabdoid_present`
		* New property: `sarcomatoid_percent`
		* New property: `sarcomatoid_present`
		* Changes made to `dysplasia_degree`
			* New permissable value: `Mild`
			* New permissable value: `Moderate`
			* New permissable value: `Severe`
		* Changes made to `dysplasia_type`
			* New permissable value: `Epithelial`
			* New permissable value: `Keratinizing`
			* New permissable value: `Nonkeratinizing`
		* Changes made to `lymph_node_involvement`
* Altered `diagnosis` Entity
	* Changes made to `deprecated`
	* Removed permissable value: `non_nodal_regional_disease`
	* Removed permissable value: `non_nodal_tumor_deposits`
	* Removed permissable value: `tumor_stage`
	* Changes made to `properties`
		* New property: `ann_arbor_b_symptoms_described`
		* Removed property: `anaplasia_present`
		* Removed property: `anaplasia_present_type`
		* Removed property: `breslow_thickness`
		* Removed property: `circumferential_resection_margin`
		* Removed property: `greatest_tumor_dimension`
		* Removed property: `gross_tumor_weight`
		* Removed property: `largest_extrapelvic_peritoneal_focus`
		* Removed property: `lymph_nodes_positive`
		* Removed property: `lymph_nodes_tested`
		* Removed property: `lymphatic_invasion_present`
		* Removed property: `lymph_node_involved_site`
		* Removed property: `non_nodal_regional_disease`
		* Removed property: `non_nodal_tumor_deposits`
		* Removed property: `percent_tumor_invasion`
		* Removed property: `perineural_invasion_present`
		* Removed property: `peripancreatic_lymph_nodes_positive`
		* Removed property: `peripancreatic_lymph_nodes_tested`
		* Removed property: `transglottic_extension`
		* Removed property: `tumor_largest_dimension_diameter`
		* Removed property: `tumor_stage`
		* Removed property: `vascular_invasion_present`
		* Removed property: `vascular_invasion_type`
		* Changes made to `ajcc_clinical_stage`
			* New permissable value: `Stage IA3`
		* Changes made to `ajcc_pathologic_m`
			* New permissable value: `M1d`
		* Changes made to `metastasis_at_diagnosis_site`
			* New permissable value: `Gastrointestinal Tract`
			* New permissable value: `Heart`
			* New permissable value: `Neck`
			* New permissable value: `Retroperitoneum`
			* New permissable value: `Urethra`
			* New permissable value: `Uterine Adnexa`
			* New permissable value: `Vertebral Canal`
			* New permissable value: `Vulva, NOS`
		* Changes made to `morphology`
			* New permissable value: `8249/6`
			* New permissable value: `8800/6`
		* Changes made to `supratentorial_localization`
			* New permissable value: `Frontal lobe`
			* New permissable value: `Occipital lobe`
			* New permissable value: `Parietal lobe`
			* New permissable value: `Temporal lobe`
* Altered `treatment` Entity
	* Changes made to `properties`
		* Changes made to `treatment_dose_units`
			* New permissable value: `mg`
		* Changes made to `treatment_type`
			* New permissable value: `Radiation, Hypofractionated`
			* New permissable value: `Radiation, Mixed Photon Beam`
			* New permissable value: `Radiation, Photon Beam`
		* Changes made to `therapeutic_agents`
			* New permissable value: `3'-dA Phosphoramidate NUC-7738`
			* New permissable value: `A2A Receptor Antagonist EOS100850`
			* New permissable value: `Adenosine A2A Receptor Antagonist CS3005`
			* New permissable value: `Adenovirus 5/F35-Human Guanylyl Cyclase C-PADRE`
			* New permissable value: `Adenovirus Serotype 26-expressing HPV16 Vaccine JNJ-63682918`
			* New permissable value: `Adenovirus Serotype 26-expressing HPV18 Vaccine JNJ-63682931`
			* New permissable value: `ALK Inhibitor TAE684`
			* New permissable value: `ALK/ROS1/Met Inhibitor TQ-B3101`
			* New permissable value: `Allogeneic Anti-BCMA CAR-transduced T-cells ALLO-715`
			* New permissable value: `Allogeneic Anti-BCMA-CAR T-cells PBCAR269A`
			* New permissable value: `Allogeneic Anti-BCMA/CS1 Bispecific CAR-T Cells`
			* New permissable value: `Allogeneic Anti-CD19 CAR T-cells ALLO-501A`
			* New permissable value: `Allogeneic Anti-CD19 Universal CAR-T Cells CTA101`
			* New permissable value: `Allogeneic Anti-CD20 CAR T-cells LUCAR-20S`
			* New permissable value: `Allogeneic Anti-CD20-CAR T-cells PBCAR20A`
			* New permissable value: `Allogeneic CD22-specific Universal CAR-expressing T-lymphocytes UCART22`
			* New permissable value: `Allogeneic CD56-positive CD3-negative Natural Killer Cells CYNK-001`
			* New permissable value: `Allogeneic CD8+ Leukemia-associated Antigens Specific T Cells NEXI-001`
			* New permissable value: `Allogeneic CRISPR-Cas9 Engineered Anti-BCMA T Cells CTX120`
			* New permissable value: `Allogeneic CRISPR-Cas9 Engineered Anti-CD70 CAR-T Cells CTX130`
			* New permissable value: `Allogeneic CS1-specific Universal CAR-expressing T-lymphocytes UCARTCS1A`
			* New permissable value: `Allogeneic Plasmacytoid Dendritic Cells Expressing Lung Tumor Antigens PDC*lung01`
			* New permissable value: `Allogeneic Third-party Suicide Gene-transduced Anti-HLA-DPB1*0401 CD4+ T-cells CTL 19`
			* New permissable value: `Allosteric ErbB Inhibitor BDTX-189`
			* New permissable value: `Alobresib`
			* New permissable value: `Alofanib`
			* New permissable value: `Alpha V Beta 8 Antagonist PF-06940434`
			* New permissable value: `Alsevalimab`
			* New permissable value: `Amivantamab`
			* New permissable value: `Andecaliximab`
			* New permissable value: `Androgen Receptor Degrader CC-94676`
			* New permissable value: `Androgen Receptor Inhibitor EPI-7386`
			* New permissable value: `Androgen Receptor/Glucocorticoid Receptor  Antagonist CB-03-10`
			* New permissable value: `Anhydrous Enol-oxaloacetate`
			* New permissable value: `Anti-5T4 Antibody-drug Conjugate ASN004`
			* New permissable value: `Anti-5T4 Antibody-drug Conjugate SYD1875`
			* New permissable value: `Anti-B7-H3/DXd Antibody-drug Conjugate DS-7300a`
			* New permissable value: `Anti-BCMA Antibody-drug Conjugate CC-99712`
			* New permissable value: `Anti-BCMA SparX Protein Plus BCMA-directed Anti-TAAG ARC T-cells CART-ddBCMA`
			* New permissable value: `Anti-BCMA/Anti-CD3 Bispecific Antibody REGN5459`
			* New permissable value: `Anti-BTLA Monoclonal Antibody TAB004`
			* New permissable value: `Anti-BTN3A Agonistic Monoclonal Antibody ICT01`
			* New permissable value: `Anti-c-Met Monoclonal Antibody HLX55`
			* New permissable value: `Anti-CCR7 Antibody-drug Conjugate JBH492`
			* New permissable value: `Anti-CD117 Monoclonal Antibody JSP191`
			* New permissable value: `Anti-CD137 Agonistic Monoclonal Antibody AGEN2373`
			* New permissable value: `Anti-CD137 Agonistic Monoclonal Antibody ATOR-1017`
			* New permissable value: `Anti-CD137 Agonistic Monoclonal Antibody LVGN6051`
			* New permissable value: `Anti-CD19 Antibody-T-cell Receptor-expressing T-cells ET019003`
			* New permissable value: `Anti-CD19 iCAR NK Cells`
			* New permissable value: `Anti-CD19/CD22 CAR NK Cells`
			* New permissable value: `Anti-CD20 Monoclonal Antibody BAT4306F`
			* New permissable value: `Anti-CD20 Monoclonal Antibody MIL62`
			* New permissable value: `Anti-CD205 Antibody-drug Conjugate OBT076`
			* New permissable value: `Anti-CD228/MMAE Antibody-drug Conjugate SGN-CD228A`
			* New permissable value: `Anti-CD25 Monoclonal Antibody RO7296682`
			* New permissable value: `Anti-CD27 Agonistic Monoclonal Antibody MK-5890`
			* New permissable value: `Anti-CD3/Anti-5T4 Bispecific Antibody GEN1044`
			* New permissable value: `Anti-CD3/Anti-GUCY2C Bispecific Antibody PF-07062119`
			* New permissable value: `Anti-CD3/CD7-Ricin Toxin A Immunotoxin`
			* New permissable value: `Anti-CD30/DM1 Antibody-drug Conjugate F0002`
			* New permissable value: `Anti-CD37 Bispecific Monoclonal Antibody GEN3009`
			* New permissable value: `Anti-CD38 Antibody-drug Conjugate STI-6129`
			* New permissable value: `Anti-CD38 Monoclonal Antibody SAR442085`
			* New permissable value: `Anti-CD38/CD28xCD3 Tri-specific Monoclonal Antibody SAR442257`
			* New permissable value: `Anti-CD39 Monoclonal Antibody SRF617`
			* New permissable value: `Anti-CD40/Anti-4-1BB Bispecific Agonist Monoclonal Antibody GEN1042`
			* New permissable value: `Anti-CD47 ADC SGN-CD47M`
			* New permissable value: `Anti-CD47 Monoclonal Antibody IMC-002`
			* New permissable value: `Anti-CD47/CD19 Bispecific Monoclonal Antibody TG-1801`
			* New permissable value: `Anti-claudin18.2 Monoclonal Antibody AB011`
			* New permissable value: `Anti-Claudin18.2 Monoclonal Antibody TST001`
			* New permissable value: `Anti-CTLA-4 Monoclonal Antibody ADG116`
			* New permissable value: `Anti-CTLA-4 Monoclonal Antibody HBM4003`
			* New permissable value: `Anti-CTLA-4 Monoclonal Antibody ONC-392`
			* New permissable value: `Anti-CTLA-4 Probody BMS-986288`
			* New permissable value: `Anti-CTLA-4/Anti-PD-1 Monoclonal Antibody Combination BCD-217`
			* New permissable value: `Anti-CTLA4 Antibody Fc Fusion Protein KN044`
			* New permissable value: `Anti-EGFR/CD16A Bispecific Antibody AFM24`
			* New permissable value: `Anti-FRA/Eribulin Antibody-drug Conjugate MORAb-202`
			* New permissable value: `Anti-GARP Monoclonal Antibody ABBV-151`
			* New permissable value: `Anti-Globo H/MMAE Antibody-drug Conjugate OBI 999`
			* New permissable value: `Anti-GPR20/DXd Antibody-drug Conjugate DS-6157a`
			* New permissable value: `Anti-gremlin-1 Monoclonal Antibody UCB6114`
			* New permissable value: `Anti-HER2 Antibody Conjugated Natural Killer Cells ACE1702`
			* New permissable value: `Anti-HER2 Antibody-drug Conjugate BAT8001`
			* New permissable value: `Anti-HER2 Antibody-drug Conjugate DP303c`
			* New permissable value: `Anti-HER2 Monoclonal Antibody B002`
			* New permissable value: `Anti-HER2 Monoclonal Antibody HLX22`
			* New permissable value: `Anti-HER2-DM1 ADC B003`
			* New permissable value: `Anti-HER2-DM1 Antibody-drug Conjugate GQ1001`
			* New permissable value: `Anti-HER2/MMAE Antibody-drug Conjugate MRG002`
			* New permissable value: `Anti-HLA-G Antibody TTX-080`
			* New permissable value: `Anti-IL-8 Monoclonal Antibody BMS-986253`
			* New permissable value: `Anti-integrin Beta-6/MMAE Antibody-drug Conjugate SGN-B6A`
			* New permissable value: `Anti-IRF4 Antisense Oligonucleotide ION251`
			* New permissable value: `Anti-LAG-3 Monoclonal Antibody IBI-110`
			* New permissable value: `Anti-latent TGF-beta 1 Monoclonal Antibody SRK-181`
			* New permissable value: `Anti-Lewis B/Lewis Y Monoclonal Antibody GNX102`
			* New permissable value: `Anti-LILRB4 Monoclonal Antibody IO-202`
			* New permissable value: `Anti-mesothelin/MMAE Antibody-drug Conjugate RC88`
			* New permissable value: `Anti-MUC16/CD3 Bispecific Antibody REGN4018`
			* New permissable value: `Anti-MUC17/CD3 BiTE Antibody AMG 199`
			* New permissable value: `Anti-NaPi2b Antibody-drug Conjugate XMT-1592`
			* New permissable value: `Anti-NRP1 Antibody ASP1948`
			* New permissable value: `Anti-OX40 Agonist Monoclonal Antibody BGB-A445`
			* New permissable value: `Anti-OX40 Hexavalent Agonist Antibody INBRX-106`
			* New permissable value: `Anti-PD-1 Antibody-interleukin-21 Mutein Fusion Protein AMG 256`
			* New permissable value: `Anti-PD-1 Monoclonal Antibody 609A`
			* New permissable value: `Anti-PD-1 Monoclonal Antibody SCT-I10A`
			* New permissable value: `Anti-PD-1/Anti-HER2 Bispecific Antibody IBI315`
			* New permissable value: `Anti-PD-1/Anti-LAG-3 Bispecific Antibody RO7247669`
			* New permissable value: `Anti-PD-1/Anti-PD-L1 Bispecific Antibody IBI318`
			* New permissable value: `Anti-PD-1/CD47 Infusion Protein HX009`
			* New permissable value: `Anti-PD-1/CTLA-4 Bispecific Antibody MEDI5752`
			* New permissable value: `Anti-PD-1/VEGF Bispecific Antibody AK112`
			* New permissable value: `Anti-PD-L1 Monoclonal Antibody IMC-001`
			* New permissable value: `Anti-PD-L1 Monoclonal Antibody RC98`
			* New permissable value: `Anti-PD-L1/Anti-4-1BB Bispecific Monoclonal Antibody GEN1046`
			* New permissable value: `Anti-PD-L1/IL-15 Fusion Protein KD033`
			* New permissable value: `Anti-PRAME T-cell Receptor/Anti-CD3 scFv Fusion Protein IMC-F106C`
			* New permissable value: `Anti-PSMA/CD3 Bispecific Antibody CCW702`
			* New permissable value: `Anti-RANKL Monoclonal Antibody GB-223`
			* New permissable value: `Anti-RANKL Monoclonal Antibody JMT103`
			* New permissable value: `Anti-Ribonucleoprotein Antibody ATRC-101`
			* New permissable value: `Anti-ROR1/PNU-159682 Derivative Antibody-drug Conjugate NBE-002`
			* New permissable value: `Anti-TIGIT Monoclonal Antibody BGB-A1217`
			* New permissable value: `Anti-TIGIT Monoclonal Antibody COM902`
			* New permissable value: `Anti-TIGIT Monoclonal Antibody SGN-TGT`
			* New permissable value: `Anti-TIM3 Monoclonal Antibody SHR-1702`
			* New permissable value: `Anti-TRAILR2/CDH17 Tetravalent Bispecific Antibody BI 905711`
			* New permissable value: `Anti-TROP2 Antibody-drug Conjugate BAT8003`
			* New permissable value: `Anti-TROP2 Antibody-drug Conjugate SKB264`
			* New permissable value: `Anti-VEGFR2 Monoclonal Antibody MSB0254`
			* New permissable value: `Antisense Oligonucleotide QR-313`
			* New permissable value: `Aprinocarsen`
			* New permissable value: `Aprutumab`
			* New permissable value: `Aryl Hydrocarbon Receptor Inhibitor IK-175`
			* New permissable value: `Aspacytarabine`
			* New permissable value: `Asunercept`
			* New permissable value: `ATR Inhibitor RP-3500`
			* New permissable value: `ATR Kinase Inhibitor M1774`
			* New permissable value: `Attenuated Measles Virus Encoding SCD Transgene TMV-018`
			* New permissable value: `Atuveciclib`
			* New permissable value: `Audencel`
			* New permissable value: `Autologous AFP Specific T Cell Receptor Transduced T Cells C-TCR055`
			* New permissable value: `Autologous Anti-BCMA CAR T-cells PHE885`
			* New permissable value: `Autologous Anti-BCMA CD8+ CAR T-cells Descartes-11`
			* New permissable value: `Autologous Anti-BCMA-CAR-4-1BB-CD3zeta-expressing T-cells C-CAR088`
			* New permissable value: `Autologous Anti-CD123 CAR-T Cells`
			* New permissable value: `Autologous Anti-CD19 CAR T-cells 19(T2)28z1xx`
			* New permissable value: `Autologous Anti-CD19 CAR-4-1BB-CD3zeta-expressing T-cells CNCT19`
			* New permissable value: `Autologous Anti-CD19 Chimeric Antigen Receptor T-cells AUTO1`
			* New permissable value: `Autologous Anti-CD19 TAC-T cells TAC01-CD19`
			* New permissable value: `Autologous Anti-CD19/CD20 Bispecific Nanobody-based CAR-T cells`
			* New permissable value: `Autologous Anti-CD19CAR-HER2t/CD22CAR-EGFRt-expressing T-cells`
			* New permissable value: `Autologous Anti-CD20 CAR Transduced CD4/CD8 Enriched T-cells MB-CART20.1`
			* New permissable value: `Autologous Anti-EGFR CAR-transduced CXCR 5-modified T-lymphocytes`
			* New permissable value: `Autologous Anti-FLT3 CAR T Cells AMG 553`
			* New permissable value: `Autologous Anti-ICAM-1-CAR-CD28-4-1BB-CD3zeta-expressing T-cells AIC100`
			* New permissable value: `Autologous Anti-kappa Light Chain CAR-CD28-expressing T-lymphocytes`
			* New permissable value: `Autologous Anti-PD-1 Antibody-activated Tumor-infiltrating Lymphocytes`
			* New permissable value: `Autologous Anti-PSMA CAR-T Cells P-PSMA-101`
			* New permissable value: `Autologous BCMA-targeted CAR T Cells CC-98633`
			* New permissable value: `Autologous Bispecific BCMA/CD19-targeted CAR-T Cells GC012F`
			* New permissable value: `Autologous Bispecific CD19/CD22-targeted CAR-T Cells GC022`
			* New permissable value: `Autologous CD19 CAR-expressing CD4+/CD8+ T-cells MB-CART19.1`
			* New permissable value: `Autologous CD19-targeted CAR T Cells CC-97540`
			* New permissable value: `Autologous CD19-targeted CAR-T Cells GC007F`
			* New permissable value: `Autologous CD19/PD-1 Bispecific CAR-T Cells`
			* New permissable value: `Autologous Clonal Neoantigen T Cells ATL001`
			* New permissable value: `Autologous CRISPR-edited Anti-CD19 CAR T Cells XYF19`
			* New permissable value: `Autologous Monocyte-derived Lysate-pulsed Dendritic Cell Vaccine PV-001-DC`
			* New permissable value: `Autologous Multi-lineage Potential Cells`
			* New permissable value: `Autologous Nectin-4/FAP-targeted CAR-T Cells`
			* New permissable value: `Autologous NKG2D CAR T-cells CYAD-02`
			* New permissable value: `Autologous Pancreatic Adenocarcinoma Lysate and mRNA-loaded Dendritic Cell Vaccine`
			* New permissable value: `Autologous Peripheral Blood Lymphocytes from Ibrutinib-treated Chronic Lymphocytic Leukemia Patients IOV-2001`
			* New permissable value: `Autologous Rapamycin-resistant Th1/Tc1 Cells RAPA-201`
			* New permissable value: `Autologous TCRm-expressing T-cells ET140203`
			* New permissable value: `Autologous Tetravalent Dendritic Cell Vaccine MIDRIX4-LUNG`
			* New permissable value: `Autologous Tumor Infiltrating Lymphocytes LN-145-S1`
			* New permissable value: `Autologous Universal CAR-expressing T-lymphocytes UniCAR02-T`
			* New permissable value: `Avdoralimab`
			* New permissable value: `Aviscumine`
			* New permissable value: `Axalimogene Filolisbac`
			* New permissable value: `Axatilimab`
			* New permissable value: `AXL Inhibitor SLC-391`
			* New permissable value: `AXL/ FLT3/VEGFR2 Inhibitor KC1036`
			* New permissable value: `Axl/Mer Inhibitor PF-07265807`
			* New permissable value: `Azintuxizumab Vedotin`
			* New permissable value: `Balstilimab`
			* New permissable value: `Bazlitoran`
			* New permissable value: `Bcl-2 Inhibitor BGB-11417`
			* New permissable value: `Bcl-2 Inhibitor LP-108`
			* New permissable value: `BCMA-CD19 Compound CAR T Cells`
			* New permissable value: `BCMA/CD3e Tri-specific T-cell Activating Construct HPN217`
			* New permissable value: `Belantamab Mafodotin`
			* New permissable value: `Belapectin`
			* New permissable value: `Belvarafenib`
			* New permissable value: `Belzutifan`
			* New permissable value: `Bersanlimab`
			* New permissable value: `Berzosertib`
			* New permissable value: `Betaglucin Gel`
			* New permissable value: `Bexmarilimab`
			* New permissable value: `Bispecific Antibody AGEN1223`
			* New permissable value: `Bispecific Antibody AMG 509`
			* New permissable value: `Bispecific Antibody GS-1423`
			* New permissable value: `BiTE Antibody AMG 910`
			* New permissable value: `Bizalimogene Ralaplasmid`
			* New permissable value: `Bomedemstat`
			* New permissable value: `BRAF Inhibitor BGB-3245`
			* New permissable value: `BRAF(V600E) Kinase Inhibitor ABM-1310`
			* New permissable value: `BTK Inhibitor HZ-A-018`
			* New permissable value: `c-Met Inhibitor ABN401`
			* New permissable value: `c-Met Inhibitor GST-HG161`
			* New permissable value: `C/EBP Beta Antagonist ST101`
			* New permissable value: `Calcium Release-activated Channel Inhibitor CM4620`
			* New permissable value: `Camidanlumab Tesirine`
			* New permissable value: `Carbon C 14-pamiparib`
			* New permissable value: `Cationic Peptide Cream Cypep-1`
			* New permissable value: `CD11b Agonist GB1275`
			* New permissable value: `CD123-CD33 Compound CAR T Cells`
			* New permissable value: `CD123-specific Targeting Module TM123`
			* New permissable value: `CD20-CD19 Compound CAR T Cells`
			* New permissable value: `CD28/ICOS Antagonist ALPN-101`
			* New permissable value: `CD44v6-specific CAR T-cells`
			* New permissable value: `CD73 Inhibitor AB680`
			* New permissable value: `CD73 Inhibitor LY3475070`
			* New permissable value: `CD80-Fc Fusion Protein ALPN-202`
			* New permissable value: `CD80-Fc Fusion Protein FPT155`
			* New permissable value: `CDK2 Inhibitor PF-07104091`
			* New permissable value: `CDK4/6 Inhibitor CS3002`
			* New permissable value: `CDK4/6 Inhibitor HS-10342`
			* New permissable value: `CDK4/6 Inhibitor TQB3616`
			* New permissable value: `CDK7 Inhibitor SY-5609`
			* New permissable value: `CDK8/19 Inhibitor SEL 120`
			* New permissable value: `Cedazuridine/Azacitidine Combination Agent ASTX030`
			* New permissable value: `Cereblon E3 Ubiquitin Ligase Modulating Agent CC-99282`
			* New permissable value: `Cetuximab Sarotalocan`
			* New permissable value: `Cevostamab`
			* New permissable value: `Chlorotoxin (EQ)-CD28-CD3zeta-CD19t-expressing CAR T-lymphocytes`
			* New permissable value: `Ciltacabtagene Autoleucel`
			* New permissable value: `Cinrebafusp Alfa`
			* New permissable value: `Cintirorgon`
			* New permissable value: `CK1alpha/CDK7/CDK9 Inhibitor BTX-A51`
			* New permissable value: `Cobolimab`
			* New permissable value: `Cofetuzumab Pelidotin`
			* New permissable value: `Coltuximab Ravtansine`
			* New permissable value: `Commensal Bacterial Strain Formulation VE800`
			* New permissable value: `Copper Cu 67 Tyr3-octreotate`
			* New permissable value: `Cord Blood Derived CAR T-Cells`
			* New permissable value: `Cosibelimab`
			* New permissable value: `Coxsackievirus V937`
			* New permissable value: `CSF1R Inhibitor ABSK021`
			* New permissable value: `CXCR4/E-selectin Antagonist GMI-1359`
			* New permissable value: `CYP11A1 Inhibitor ODM-209`
			* New permissable value: `CYP17/CYP11B2 Inhibitor LAE001`
			* New permissable value: `Daratumumab and Hyaluronidase-fihj`
			* New permissable value: `Decitabine and Cedazuridine`
			* New permissable value: `Delolimogene Mupadenorepvec`
			* New permissable value: `Dendrimer-conjugated Bcl-2/Bcl-XL Inhibitor AZD0466`
			* New permissable value: `Dengue Virus Adjuvant PV-001-DV`
			* New permissable value: `Dilpacimab`
			* New permissable value: `DNA-PK inhibitor AZD7648`
			* New permissable value: `DNA-PK/PI3K-delta Inhibitor BR101801`
			* New permissable value: `DNMT1 Inhibitor NTX-301`
			* New permissable value: `Dociparstat sodium`
			* New permissable value: `Dostarlimab`
			* New permissable value: `Doxorubicin Prodrug/Prodrug-activating Biomaterial SQ3370`
			* New permissable value: `DTRMWXHS-12/Everolimus/Pomalidomide Combination Agent DTRM-555`
			* New permissable value: `Edicotinib`
			* New permissable value: `Edodekin alfa`
			* New permissable value: `Eftozanermin Alfa`
			* New permissable value: `EGFR Inhibitor TY-9591`
			* New permissable value: `EGFR Mutant-selective  Inhibitor TQB3804`
			* New permissable value: `EGFR/EGFRvIII Inhibitor WSD0922-FU`
			* New permissable value: `EGFR/HER2 Inhibitor DZD9008`
			* New permissable value: `EGFR/TGFb Fusion Monoclonal Antibody BCA101`
			* New permissable value: `EGFR/VEGFR/RET Inhibitor HA121-28`
			* New permissable value: `Emibetuzumab`
			* New permissable value: `Enadenotucirev-expressing FAP/CD3 Bispecific FAP-TAc NG-641`
			* New permissable value: `Encapsulated Rapamycin`
			* New permissable value: `Encelimab`
			* New permissable value: `Endothelin B Receptor Blocker ENB 003`
			* New permissable value: `Engineered Red Blood Cells Co-expressing 4-1BBL and IL-15TP RTX-240`
			* New permissable value: `Engineered Toxin Body Targeting CD38 TAK-169`
			* New permissable value: `EP2/EP4 Antagonist TPST-1495`
			* New permissable value: `EP4 Antagonist INV-1120`
			* New permissable value: `Epcoritamab`
			* New permissable value: `EphA2-targeting Bicycle Toxin Conjugate BT5528`
			* New permissable value: `ERK1/2 Inhibitor HH2710`
			* New permissable value: `ERK1/2 Inhibitor JSI-1187`
			* New permissable value: `Etigilimab`
			* New permissable value: `Exicorilant`
			* New permissable value: `Extended Release Metformin Hydrochloride`
			* New permissable value: `Ezabenlimab`
			* New permissable value: `EZH1/2 Inhibitor HH2853`
			* New permissable value: `EZH2 inhibitor CPI-0209`
			* New permissable value: `Fadraciclib`
			* New permissable value: `FAK/ALK/ROS1 Inhibitor APG-2449`
			* New permissable value: `FAP/4-1BB-targeting DARPin MP0310`
			* New permissable value: `FAP/4-1BB-targeting Fusion Protein RO7122290`
			* New permissable value: `Fas Ligand-treated Allogeneic Mobilized Peripheral Blood Cells`
			* New permissable value: `Favezelimab`
			* New permissable value: `Fc-engineered Anti-CD40 Agonist Antibody 2141-V11`
			* New permissable value: `Feladilimab`
			* New permissable value: `Felzartamab`
			* New permissable value: `Fenretinide Phospholipid Suspension ST-001`
			* New permissable value: `FGFR Inhibitor CPL304110`
			* New permissable value: `FGFR/CSF-1R Inhibitor 3D185`
			* New permissable value: `FGFR2 Inhibitor RLY-4008`
			* New permissable value: `Fianlimab`
			* New permissable value: `Fimaporfin A`
			* New permissable value: `Flt3 Ligand/Anti-CTLA-4 Antibody/IL-12 Engineered Oncolytic Vaccinia Virus RIVAL-01`
			* New permissable value: `FLT3/FGFR Dual Kinase Inhibitor MAX-40279`
			* New permissable value: `FLT3/KIT/CSF1R Inhibitor NMS-03592088`
			* New permissable value: `Fluorine F 18 Ara-G`
			* New permissable value: `Foritinib Succinate`
			* New permissable value: `Fosgemcitabine Palabenamide`
			* New permissable value: `Fosifloxuridine Nafalbenamide`
			* New permissable value: `Futibatinib`
			* New permissable value: `G Protein-coupled Estrogen Receptor Agonist LNS8801`
			* New permissable value: `Gallium-based Bone Resorption Inhibitor AP-002`
			* New permissable value: `GBM Antigens and Alloantigens Immunotherapeutic Vaccine`
			* New permissable value: `Genetically Modified Interleukin-12 Transgene-encoding Bifidobacterium longum`
			* New permissable value: `Giloralimab`
			* New permissable value: `Giredestrant`
			* New permissable value: `Glofitamab`
			* New permissable value: `Glutaminase Inhibitor IPN60090`
			* New permissable value: `Glutamine Antagonist DRP-104`
			* New permissable value: `HER2 Inhibitor DZD1516`
			* New permissable value: `HER2 Tri-specific Natural Killer Cell Engager DF1001`
			* New permissable value: `HER2-directed TLR8 Agonist SBT6050`
			* New permissable value: `HIF2a RNAi ARO-HIF2`
			* New permissable value: `HPV 16 E6/E7-encoding Arenavirus Vaccine HB-202`
			* New permissable value: `HPV E6/E7-encoding Arenavirus Vaccine HB-201`
			* New permissable value: `HPV6/11-targeted DNA Plasmid Vaccine INO-3107`
			* New permissable value: `Hsp90 Inhibitor TQB3474`
			* New permissable value: `Hsp90-targeted Photosensitizer HS-201`
			* New permissable value: `Hyaluronidase-zzxf/Pertuzumab/Trastuzumab`
			* New permissable value: `Iadademstat`
			* New permissable value: `Ianalumab`
			* New permissable value: `Idetrexed`
			* New permissable value: `IDH1 Mutant Inhibitor LY3410738`
			* New permissable value: `IDO/TDO Inhibitor LY-01013`
			* New permissable value: `IDO1/TDO2 Inhibitor M4112`
			* New permissable value: `Ieramilimab`
			* New permissable value: `IL-12sc, IL-15sushi, IFNa and GM-CSF mRNA-based Immunotherapeutic Agent SAR441000`
			* New permissable value: `Ilginatinib`
			* New permissable value: `Imaradenant`
			* New permissable value: `Imgatuzumab`
			* New permissable value: `Imifoplatin`
			* New permissable value: `Inactivated Oncolytic Virus Particle GEN0101`
			* New permissable value: `Indatuximab Ravtansine`
			* New permissable value: `Individualized MVA-based Vaccine TG4050`
			* New permissable value: `Indusatumab Vedotin`
			* New permissable value: `Inebilizumab`
			* New permissable value: `iNKT Cell Agonist ABX196`
			* New permissable value: `Interleukin-12-Fc Fusion Protein DF6002`
			* New permissable value: `Interleukin-15 Agonist Fusion Protein SHR1501`
			* New permissable value: `Interleukin-15 Fusion Protein BJ-001`
			* New permissable value: `Interleukin-15/Interleukin-15 Receptor Alpha Complex-Fc Fusion Protein XmAb24306`
			* New permissable value: `Interleukin-15/Interleukin-15 Receptor Alpha Sushi+ Domain Fusion Protein SO-C101`
			* New permissable value: `Iodine I 131 Apamistamab`
			* New permissable value: `Iodine I 131 IPA`
			* New permissable value: `iPSC-derived CD16/IL-15RF-expressing Anti-CD19 CAR-NK Cells FT596`
			* New permissable value: `Irinotecan Sucrosofate`
			* New permissable value: `Iroplact`
			* New permissable value: `Irradiated Allogeneic Human Lung Cancer Cells Expressing OX40L-Ig Vaccine HS-130`
			* New permissable value: `Istiratumab`
			* New permissable value: `Ivaltinostat`
			* New permissable value: `Ivuxolimab`
			* New permissable value: `JAK Inhibitor`
			* New permissable value: `Kanitinib`
			* New permissable value: `KRAS G12C Inhibitor GDC-6036`
			* New permissable value: `KRAS G12C Inhibitor LY3499446`
			* New permissable value: `KRASG12C Inhibitor JNJ-74699157`
			* New permissable value: `Ladiratuzumab Vedotin`
			* New permissable value: `LAIR-2 Fusion Protein NC410`
			* New permissable value: `Landogrozumab`
			* New permissable value: `Laprituximab Emtansine`
			* New permissable value: `Larotinib Mesylate`
			* New permissable value: `Lerociclib`
			* New permissable value: `Letetresgene Autoleucel`
			* New permissable value: `Letolizumab`
			* New permissable value: `Lifileucel`
			* New permissable value: `Lifirafenib`
			* New permissable value: `Lilotomab`
			* New permissable value: `Linperlisib`
			* New permissable value: `Lipid Nanoparticle Encapsulating Glutathione S-transferase P siRNA NBF-006`
			* New permissable value: `Liposomal Bcl-2 Antisense Oligonucleotide  BP1002`
			* New permissable value: `Liposome-encapsulated TAAs mRNA Vaccine W_ova1`
			* New permissable value: `LMP2-specific T Cell Receptor-transduced Autologous T-lymphocytes`
			* New permissable value: `LMP7 Inhibitor M3258`
			* New permissable value: `Lodapolimab`
			* New permissable value: `LRP5 Antagonist BI 905681`
			* New permissable value: `LSD1 Inhibitor SYHA1807`
			* New permissable value: `Lumretuzumab`
			* New permissable value: `Lupartumab Amadotin`
			* New permissable value: `Lutetium Lu 177-DTPA-omburtamab`
			* New permissable value: `Maackia amurensis Seed Lectin`
			* New permissable value: `Macrocycle-bridged STING Agonist E7766`
			* New permissable value: `MAGE-A1-specific T Cell Receptor-transduced Autologous T-cells`
			* New permissable value: `Magrolimab`
			* New permissable value: `Manelimab`
			* New permissable value: `MCL-1 Inhibitor ABBV-467`
			* New permissable value: `MDM2 Inhibitor AMGMDS3`
			* New permissable value: `MEK 1/2 Inhibitor FCN-159`
			* New permissable value: `MEK Inhibitor HL-085`
			* New permissable value: `Menin-MLL Interaction Inhibitor SNDX-5613`
			* New permissable value: `MET x MET Bispecific Antibody REGN5093`
			* New permissable value: `Metarrestin`
			* New permissable value: `Methylcantharidimide`
			* New permissable value: `Mevociclib`
			* New permissable value: `Mezagitamab`
			* New permissable value: `Microbiome GEN-001`
			* New permissable value: `Microbiome-derived Peptide Vaccine EO2401`
			* New permissable value: `Milataxel`
			* New permissable value: `Miptenalimab`
			* New permissable value: `Miransertib`
			* New permissable value: `Mirdametinib`
			* New permissable value: `Mirzotamab Clezutoclax`
			* New permissable value: `Mitazalimab`
			* New permissable value: `Mobocertinib`
			* New permissable value: `Modakafusp Alfa`
			* New permissable value: `Modified Vaccinia Ankara-vectored HPV16/18 Vaccine JNJ-65195208`
			* New permissable value: `Motixafortide`
			* New permissable value: `MUC-1/WT1 Peptide-primed Autologous Dendritic Cells`
			* New permissable value: `Multi-epitope HER2 Peptide Vaccine TPIV100`
			* New permissable value: `Murizatoclax`
			* New permissable value: `Muscadine Grape Extract`
			* New permissable value: `MVA-BN Smallpox Vaccine`
			* New permissable value: `N-dihydrogalactochitosan`
			* New permissable value: `Nagrestipen`
			* New permissable value: `Naratuximab Emtansine`
			* New permissable value: `Navicixizumab`
			* New permissable value: `Nogapendekin Alfa`
			* New permissable value: `Numidargistat`
			* New permissable value: `Nurulimab`
			* New permissable value: `Odronextamab`
			* New permissable value: `Ofranergene Obadenovec`
			* New permissable value: `Oligo-fucoidan`
			* New permissable value: `Olinvacimab`
			* New permissable value: `Olvimulogene Nanivacirepvec`
			* New permissable value: `Onatasertib`
			* New permissable value: `Oncolytic Adenovirus ORCA-010`
			* New permissable value: `Oncolytic Herpes Simplex Virus-1 ONCR-177`
			* New permissable value: `Oncolytic HSV-1 Expressing IL-12 and Anti-PD-1 Antibody T3011`
			* New permissable value: `Oncolytic Measles Virus Encoding Helicobacter pylori Neutrophil-activating Protein`
			* New permissable value: `Ontorpacept`
			* New permissable value: `Onvatilimab`
			* New permissable value: `Opolimogene Capmilisbac`
			* New permissable value: `Opucolimab`
			* New permissable value: `Orelabrutinib`
			* New permissable value: `Orvacabtagene Autoleucel`
			* New permissable value: `Oxaliplatin Eluting Beads`
			* New permissable value: `p97 Inhibitor CB-5339`
			* New permissable value: `p97 Inhibitor CB-5339 Tosylate`
			* New permissable value: `Pacmilimab`
			* New permissable value: `Pamrevlumab`
			* New permissable value: `Pan-KRAS Inhibitor BI 1701963`
			* New permissable value: `Pan-mutation-selective EGFR Inhibitor CLN-081`
			* New permissable value: `Pan-TRK Inhibitor NOV1601`
			* New permissable value: `Panulisib`
			* New permissable value: `PARP 1/2 Inhibitor IMP4297`
			* New permissable value: `PARP Inhibitor NMS-03305293`
			* New permissable value: `PARP7 Inhibitor RBN-2397`
			* New permissable value: `Parsaclisib`
			* New permissable value: `Parsaclisib Hydrochloride`
			* New permissable value: `Partially Engineered T-regulatory Cell Donor Graft TRGFT-201`
			* New permissable value: `Patritumab Deruxtecan`
			* New permissable value: `Paxalisib`
			* New permissable value: `PD-L1 Inhibitor GS-4224`
			* New permissable value: `PD-L1/4-1BB/HSA Trispecific Fusion Protein NM21-1480`
			* New permissable value: `Pegvorhyaluronidase Alfa`
			* New permissable value: `Pegylated SN-38 Conjugate PLX038`
			* New permissable value: `Pelabresib`
			* New permissable value: `Peposertib`
			* New permissable value: `Personalized and Adjusted Neoantigen Peptide Vaccine PANDA-VAC`
			* New permissable value: `Personalized Neoantigen DNA Vaccine GNOS-PV01`
			* New permissable value: `Personalized Neoantigen DNA Vaccine GNOS-PVO2`
			* New permissable value: `Photodynamic Compound TLD-1433`
			* New permissable value: `Pimitespib`
			* New permissable value: `Pimurutamab`
			* New permissable value: `Pinatuzumab Vedotin`
			* New permissable value: `Pixatimod`
			* New permissable value: `Plamotamab`
			* New permissable value: `Plasmid DNA Vaccine pING-hHER3FL`
			* New permissable value: `pNGVL4a-CRT-E6E7L2 DNA Vaccine`
			* New permissable value: `pNGVL4a-Sig/E7(detox)/HSP70 DNA and HPV16 L2/E6/E7 Fusion Protein TA-CIN Vaccine PVX-2`
			* New permissable value: `Polymer-conjugated IL-15 Receptor Agonist NKTR-255`
			* New permissable value: `Porcupine Inhibitor XNW7201`
			* New permissable value: `PPAR Alpha Antagonist TPST-1120`
			* New permissable value: `Pralsetinib`
			* New permissable value: `Praluzatamab Ravtansine`
			* New permissable value: `PRMT5 Inhibitor PRT811`
			* New permissable value: `Prolgolimab`
			* New permissable value: `Prostaglandin E2 EP4 Receptor Inhibitor AN0025`
			* New permissable value: `Protein Tyrosine Kinase 2 Inhibitor IN10018`
			* New permissable value: `Pyruvate Kinase M2 Isoform Activator TP-1454`
			* New permissable value: `Racemetyrosine/Methoxsalen/Phenytoin/Sirolimus SM-88`
			* New permissable value: `Radgocitabine`
			* New permissable value: `Radgocitabine Hydrochloride`
			* New permissable value: `Ragifilimab`
			* New permissable value: `Recombinant Bacterial Minicells VAX014`
			* New permissable value: `Recombinant Erwinia asparaginase JZP-458`
			* New permissable value: `Recombinant Human Papillomavirus 11-valent Vaccine`
			* New permissable value: `Recombinant Human TRAIL-Trimer Fusion Protein SCB-313`
			* New permissable value: `Recombinant Humanized Anti-HER-2 Bispecific Monoclonal Antibody MBS301`
			* New permissable value: `Redaporfin`
			* New permissable value: `RET/SRC Inhibitor TPX-0046`
			* New permissable value: `Retifanlimab`
			* New permissable value: `Revdofilimab`
			* New permissable value: `Rezivertinib`
			* New permissable value: `Ripertamab`
			* New permissable value: `Roblitinib`
			* New permissable value: `ROBO1-targeted BiCAR-NKT Cells`
			* New permissable value: `Rocakinogene Sifuplasmid`
			* New permissable value: `Roducitabine`
			* New permissable value: `Rolinsatamab Talirine`
			* New permissable value: `Roneparstat`
			* New permissable value: `Ropeginterferon Alfa-2B`
			* New permissable value: `Ropocamptide`
			* New permissable value: `Rosopatamab`
			* New permissable value: `RSK1-4 Inhibitor PMD-026`
			* New permissable value: `Ruthenium-based Small Molecule Therapeutic BOLD-100`
			* New permissable value: `Ruxotemitide`
			* New permissable value: `Sabatolimab`
			* New permissable value: `Samrotamab Vedotin`
			* New permissable value: `Samuraciclib`
			* New permissable value: `Sasanlimab`
			* New permissable value: `SDF-1 Receptor Antagonist PTX-9908`
			* New permissable value: `Selective Estrogen Receptor Degrader LX-039`
			* New permissable value: `Selective Estrogen Receptor Degrader LY3484356`
			* New permissable value: `Serclutamab Talirine`
			* New permissable value: `SERD ZN-c5`
			* New permissable value: `Serplulimab`
			* New permissable value: `Shenqi Fuzheng Injection SQ001`
			* New permissable value: `SHP2 Inhibitor RLY-1971`
			* New permissable value: `Simlukafusp Alfa`
			* New permissable value: `Simmitinib`
			* New permissable value: `Simurosertib`
			* New permissable value: `Siremadlin`
			* New permissable value: `SIRPa-4-1BBL Fusion Protein DSP107`
			* New permissable value: `SIRPa-Fc-CD40L Fusion Protein SL-172154`
			* New permissable value: `Sotigalimab`
			* New permissable value: `Sotorasib`
			* New permissable value: `Spanlecortemlocel`
			* New permissable value: `SRPK1/ABCG2 Inhibitor SCO-101`
			* New permissable value: `STING Agonist BMS-986301`
			* New permissable value: `STING Agonist GSK3745417`
			* New permissable value: `STING Agonist IMSA101`
			* New permissable value: `STING Agonist SB 11285`
			* New permissable value: `STING Agonist TAK-676`
			* New permissable value: `STING-expressing E. coli SYNB1891`
			* New permissable value: `Sugemalimab`
			* New permissable value: `Superoxide Dismutase Mimetic GC4711`
			* New permissable value: `Synthetic Plumbagin PCUR-101`
			* New permissable value: `Tafasitamab`
			* New permissable value: `Taletrectinib`
			* New permissable value: `Tamrintamab Pamozirine`
			* New permissable value: `Tankyrase Inhibitor STP1002`
			* New permissable value: `Tapotoclax`
			* New permissable value: `Tasadenoturev`
			* New permissable value: `Tebentafusp`
			* New permissable value: `Teclistamab`
			* New permissable value: `Tefinostat`
			* New permissable value: `Telaglenastat`
			* New permissable value: `Telaglenastat Hydrochloride`
			* New permissable value: `Telisotuzumab`
			* New permissable value: `Tepoditamab`
			* New permissable value: `TGF-beta Receptor 1 Kinase Inhibitor SH3051`
			* New permissable value: `TGF-beta Receptor 1 Kinase Inhibitor YL-13027`
			* New permissable value: `Therapeutic Cancer Vaccine ATP128`
			* New permissable value: `Thorium Th 227 Anetumab Corixetan`
			* New permissable value: `Thorium Th 227 Anti-HER2 Monoclonal Antibody BAY2701439`
			* New permissable value: `Thorium Th 227 Anti-PSMA Monoclonal Antibody BAY 2315497`
			* New permissable value: `Thymidylate Synthase Inhibitor CX1106`
			* New permissable value: `TIGIT Inhibitor M6223`
			* New permissable value: `Tilogotamab`
			* New permissable value: `Tiomolibdate Choline`
			* New permissable value: `Tirbanibulin`
			* New permissable value: `TLR7 agonist BNT411`
			* New permissable value: `TLR7 Agonist LHC165`
			* New permissable value: `TM4SF1-CAR/EpCAM-CAR-expressing Autologous T Cells`
			* New permissable value: `Tolebrutinib`
			* New permissable value: `Topotecan Sustained-release Episcleral Plaque`
			* New permissable value: `Trastuzumab Deruxtecan`
			* New permissable value: `Trastuzumab Monomethyl Auristatin F`
			* New permissable value: `Trastuzumab-TLR 7/8 Agonist BDC-1001`
			* New permissable value: `Tris-acryl Gelatin Microspheres`
			* New permissable value: `TRK Inhibitor TQB3558`
			* New permissable value: `Troriluzole`
			* New permissable value: `Tyrosine Kinase Inhibitor TL-895`
			* New permissable value: `Upifitamab`
			* New permissable value: `Urabrelimab`
			* New permissable value: `Ursolic Acid`
			* New permissable value: `Uzansertib`
			* New permissable value: `Valecobulin`
			* New permissable value: `Valemetostat`
			* New permissable value: `Vesencumab`
			* New permissable value: `Vibecotamab`
			* New permissable value: `Vibostolimab`
			* New permissable value: `Vorasidenib`
			* New permissable value: `Vosilasarm`
			* New permissable value: `Vulinacimab`
			* New permissable value: `Wee1 Inhibitor ZN-c3`
			* New permissable value: `Wee1 Kinase Inhibitor Debio 0123`
			* New permissable value: `Xevinapant`
			* New permissable value: `Xiliertinib`
			* New permissable value: `Xisomab 3G3`
			* New permissable value: `Yttrium Y 90 Tabituximab Barzuxetan`
			* New permissable value: `Zandelisib`
			* New permissable value: `Zanidatamab`
			* New permissable value: `Zelavespib`
			* New permissable value: `Zorifertinib`
			* New permissable value: `Zotatifin`
			* New permissable value: `Zotiraciclib Citrate`
			* Removed permissable value: `5-Fluorouracil (5-FU)`
			* Removed permissable value: `5-Fluorouracil (5-FU) ; Leucovorin ; Oxaliplatin (Eloxatin)`
			* Removed permissable value: `5-Fluorouracil (5-FU) with the vitamin-like drug leucovorin (also called folinic acid) or levo-leucovorin`
			* Removed permissable value: `5-Fluorouracil (5-FU), Leucovorin, Oxaliplatin (Eloxatin)`
			* Removed permissable value: `Abraxane/Avastin`
			* Removed permissable value: `AC>paclitaxel`
			* Removed permissable value: `Afatinib (BIBW 2992)`
			* Removed permissable value: `Anastrozole (Arimidex)`
			* Removed permissable value: `Avastin/Xeloda`
			* Removed permissable value: `AZD2171 (cediranib)`
			* Removed permissable value: `BCNU`
			* Removed permissable value: `Belinostat (PXD-101)`
			* Removed permissable value: `Bevacizumab (Avastin)`
			* Removed permissable value: `Bevacizumab (rhuMAb VEGF)`
			* Removed permissable value: `Bicalutamide (Casodex)`
			* Removed permissable value: `BMS-936558 (Nivolumab, MDX-1106)`
			* Removed permissable value: `Cabazitaxel (Jevtana)`
			* Removed permissable value: `Capecitabine (Xeloda)`
			* Removed permissable value: `Carboplatin, etoposide, dexamethason`
			* Removed permissable value: `Carboplatin, paclitaxel, nivolumab`
			* Removed permissable value: `Carboplatin, paclitaxel, vinorelbine, nivolumab`
			* Removed permissable value: `Carboplatin,Cyclophosphamide,Doxorubicin (Adriamycin),Etoposide (VP-16)`
			* Removed permissable value: `Carboplatin/Abraxane`
			* Removed permissable value: `Carmustine (BCNU)`
			* Removed permissable value: `Cetuximab (Erbitux)`
			* Removed permissable value: `Cisplatin and 5-Fluorouracil`
			* Removed permissable value: `Cisplatin, etoposide`
			* Removed permissable value: `Cisplatin, etoposide, carboplatin, taxol, nivolumab); surgery`
			* Removed permissable value: `Cyclophosphamide (Cytoxan)`
			* Removed permissable value: `Cyclophosphamide,Vincristine; irinotecan; temozolomide (VIT)`
			* Removed permissable value: `DA-EPOCH-R`
			* Removed permissable value: `DD4A | chemotherapy vincristine, actinomycin, doxorubicin, cyclophosphamide, and etoposide`
			* Removed permissable value: `Docetaxel (Taxotere)`
			* Removed permissable value: `Doxorubicin (Adriamycin)`
			* Removed permissable value: `Doxorubicin (Adriamycin),Etoposide (VP-16),Vincristine,Vincristine; actinomycin-D; cyclophosphamide (VAC)`
			* Removed permissable value: `Doxorubicin, liposome encapsulated (Liposome Company)`
			* Removed permissable value: `E7389 (Eribulin; Halichondrin B Analog)`
			* Removed permissable value: `EMD 121974 (Cilengitide)`
			* Removed permissable value: `Epirubicin (Ellence)`
			* Removed permissable value: `Estradiol mustard`
			* Removed permissable value: `Etoposide (VP-16)`
			* Removed permissable value: `Everolimus (RAD-001)`
			* Removed permissable value: `Exemestane (Aromasin)`
			* Removed permissable value: `FOLFIRI: 5-FU, Irinotecan, Leucovorin`
			* Removed permissable value: `Folfirinox/ mFolfirinox`
			* Removed permissable value: `FOLFOX`
			* Removed permissable value: `Fulvestrant (Faslodex)`
			* Removed permissable value: `GDC-0449 (Vismodegib)`
			* Removed permissable value: `Gemcitabine (Gemzar)`
			* Removed permissable value: `Gemcitabine abraxane (nab-paclitaxel)`
			* Removed permissable value: `High-dose methotrexate; doxorubicin; cisplatin (MAP)`
			* Removed permissable value: `Ipilimumab (BMS-734016; MDX-010 Transfectoma-derived)`
			* Removed permissable value: `Irinotecan (CPT-11, Camptosar)`
			* Removed permissable value: `Ixabepilone (BMS 247550, Ixempra)`
			* Removed permissable value: `Letrozole (Femara)`
			* Removed permissable value: `Liposomal doxorubicin`
			* Removed permissable value: `Liposomal Doxorubicin (Doxil)`
			* Removed permissable value: `MDV3100 (Enzalutamide)`
			* Removed permissable value: `Megestrol acetate (Megace)`
			* Removed permissable value: `Methoxyamine hydrochloride (TRC102)`
			* Removed permissable value: `MK-2206`
			* Removed permissable value: `None`
			* Removed permissable value: `OSI-774 (erlotinib; Tarceva)`
			* Removed permissable value: `OXALIplatin (Eloxatin)`
			* Removed permissable value: `Oxaliplatin (Eloxatin) and 5-Fluorouracil (5-FU)`
			* Removed permissable value: `Paclitaxel (Taxol)`
			* Removed permissable value: `Paclitaxel protein-bound particles (albumin-bound)`
			* Removed permissable value: `PD-0332991`
			* Removed permissable value: `Pembrolizumab (Keytruda)`
			* Removed permissable value: `Pemetrexed (Alimta; LY231514)`
			* Removed permissable value: `PF-03084014`
			* Removed permissable value: `PS-341 (bortezomib; Velcade)`
			* Removed permissable value: `R-CHOP`
			* Removed permissable value: `Ramucirumab (Cyramza)`
			* Removed permissable value: `Ramucirumab (IMC-1121B)`
			* Removed permissable value: `Regorafenib (Stivarga)`
			* Removed permissable value: `Sorafenib (BAY 43-9006; Nexavar)`
			* Removed permissable value: `Sorafenib (Nexavar)`
			* Removed permissable value: `STI571 (imatinib, Gleevec)`
			* Removed permissable value: `Sunitinib malate (SU011248 L-malate; Sutent)`
			* Removed permissable value: `Tamoxifen (Nolvadex)`
			* Removed permissable value: `Taxol (Old NSC)`
			* Removed permissable value: `Temsirolimus (CCI-779)`
			* Removed permissable value: `Trastuzumab (Herceptin)`
			* Removed permissable value: `Vincristin; Dactinomycin; Cyclophosphamide (VAC)`
			* Removed permissable value: `Vincristine, actinomycin-D, cyclophosphamide (VAC); Vincristine, irinotecan, temozolomide (VIT); Etoposide; Ifosfamide`
			* Removed permissable value: `Vincristine, doxorubicin, cyclophosphamide, ifosfamide, etoposide (VDC/IE)`
			* Removed permissable value: `Vincristine; actinomycin-D; cyclophosphamide; vincristine; irinotecan (VAC/VI)`
			* Removed permissable value: `Vincristine; doxorubicin; cyclophosphamide; ifosfamide; etoposide (VDC/IE)`
			* Removed permissable value: `Vincristine; irinotecan; temozolomide (VIT)`
			* Removed permissable value: `XL184 (Cabozantinib s-malate)`
			* Removed permissable value: `Zoledronic acid (zolendronate, Zometa)`
			* Removed property: `deprecated_enum`
* Altered `case` Entity
	* Changes made to `properties`
		* Changes made to `disease_type`
			* Removed permissable value: `Acute Lymphoblastic Leukemia`
			* Removed permissable value: `Acute Myeloid Leukemia`
			* Removed permissable value: `Adrenocortical Carcinoma`
			* Removed permissable value: `Bladder Urothelial Carcinoma`
			* Removed permissable value: `Brain Lower Grade Glioma`
			* Removed permissable value: `Breast Invasive Carcinoma`
			* Removed permissable value: `Burkitt Lymphoma`
			* Removed permissable value: `Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma`
			* Removed permissable value: `Cholangiocarcinoma`
			* Removed permissable value: `Chronic Lymphocytic Leukemia`
			* Removed permissable value: `Clear Cell Sarcoma of the Kidney`
			* Removed permissable value: `Colon Adenocarcinoma`
			* Removed permissable value: `Esophageal Carcinoma`
			* Removed permissable value: `Glioblastoma Multiforme`
			* Removed permissable value: `Head and Neck Squamous Cell Carcinoma`
			* Removed permissable value: `High-Risk Wilms Tumor`
			* Removed permissable value: `HIV+ Tumor Molecular Characterization Project - Cervical Cancer`
			* Removed permissable value: `HIV+ Tumor Molecular Characterization Project - Lung Cancer`
			* Removed permissable value: `Kidney Chromophobe`
			* Removed permissable value: `Kidney Renal Clear Cell Carcinoma`
			* Removed permissable value: `Kidney Renal Papillary Cell Carcinoma`
			* Removed permissable value: `Liver Hepatocellular Carcinoma`
			* Removed permissable value: `Lung Adenocarcinoma`
			* Removed permissable value: `Lung Squamous Cell Carcinoma`
			* Removed permissable value: `Lymphoid Neoplasm Diffuse Large B-cell Lymphoma`
			* Removed permissable value: `Mesothelioma`
			* Removed permissable value: `Multiple Myeloma`
			* Removed permissable value: `Neuroblastoma`
			* Removed permissable value: `Osteosarcoma`
			* Removed permissable value: `Ovarian Serous Cystadenocarcinoma`
			* Removed permissable value: `Pancreatic Adenocarcinoma`
			* Removed permissable value: `Pheochromocytoma and Paraganglioma`
			* Removed permissable value: `Prostate Adenocarcinoma`
			* Removed permissable value: `Rectum Adenocarcinoma`
			* Removed permissable value: `Rhabdoid Tumor`
			* Removed permissable value: `Sarcoma`
			* Removed permissable value: `Skin Cutaneous Melanoma`
			* Removed permissable value: `Stomach Adenocarcinoma`
			* Removed permissable value: `Testicular Germ Cell Tumors`
			* Removed permissable value: `Thymoma`
			* Removed permissable value: `Thyroid Carcinoma`
			* Removed permissable value: `Uterine Carcinosarcoma`
			* Removed permissable value: `Uterine Corpus Endometrial Carcinoma`
			* Removed permissable value: `Uveal Melanoma`
			* Removed property: `deprecated_enum`
		* Changes made to `primary_site`
			* Removed permissable value: `Adrenal Gland`
			* Removed permissable value: `Bile Duct`
			* Removed permissable value: `Blood`
			* Removed permissable value: `Bone`
			* Removed permissable value: `Bone Marrow`
			* Removed permissable value: `Cervix`
			* Removed permissable value: `Colorectal`
			* Removed permissable value: `Eye`
			* Removed permissable value: `Head and Neck`
			* Removed permissable value: `Liver`
			* Removed permissable value: `Lung`
			* Removed permissable value: `Lymph Nodes`
			* Removed permissable value: `Nervous System`
			* Removed permissable value: `Not Applicable`
			* Removed permissable value: `Pleura`
			* Removed permissable value: `Prostate`
			* Removed permissable value: `Soft Tissue`
			* Removed permissable value: `Thyroid`
			* Removed permissable value: `Uterus`
			* Removed property: `deprecated_enum`
* Altered `rna_expression_workflow` Entity
	* Changes made to `properties`
		* Changes made to `workflow_type`
			* New permissable value: `STAR - Smart-Seq2 Raw Counts`
			* New permissable value: `STAR - Smart-Seq2 Filtered Counts`
			* Removed permissable value: `STAR - Smart-Seq2 Gene Counts`
			* Removed permissable value: `STAR - Smart-Seq2 GeneFull Counts`
* Altered `molecular_test` Entity
	* Changes made to `properties`
		* New property: `slides`
		* Changes made to `antigen`
			* New permissable value: `Ki67`
		* Changes made to `gene_symbol`
			* New permissable value: `CHGA`
			* New permissable value: `SYP`
		* Changes made to `molecular_consequence`
			* New permissable value: `Exon Variant`
		* Changes made to `second_gene_symbol`
			* New permissable value: `CHGA`
			* New permissable value: `SYP`
* Altered `aligned_reads` Entity
	* Changes made to `properties`
		* Removed property: `proportion_coverage_10X`
		* Removed property: `proportion_coverage_30X`
		* Changes made to `experimental_strategy`
			* New permissable value: `HiChIP`
			* New permissable value: `scATAC-Seq`

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
