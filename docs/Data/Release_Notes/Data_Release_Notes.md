# Data Release Notes

| Version | Date |
|---|---|
| [v38.0](Data_Release_Notes.md#data-release-380) | August XX, 2023 |
| [v37.0](Data_Release_Notes.md#data-release-370) | March 29, 2023 |
| [v36.0](Data_Release_Notes.md#data-release-360) | December 12, 2022 |
| [v35.0](Data_Release_Notes.md#data-release-350) | September 28, 2022 |
| [v34.0](Data_Release_Notes.md#data-release-340) | July 27, 2022 |
| [v33.1](Data_Release_Notes.md#data-release-331) | May 31, 2022 |
| [v33.0](Data_Release_Notes.md#data-release-330) | May 3, 2022 |
| [v32.0](Data_Release_Notes.md#data-release-320) | March 29, 2022 |
| [v31.0](Data_Release_Notes.md#data-release-310) | October 29, 2021 |
| [v30.0](Data_Release_Notes.md#data-release-300) | September 23, 2021 |
| [v29.0](Data_Release_Notes.md#data-release-290) | March 31, 2021 |
| [v28.0](Data_Release_Notes.md#data-release-280) | February 2, 2021 |
| [v27.0-fix](Data_Release_Notes.md#data-release-270-bug-fix) | November 9, 2020 |
| [v27.0](Data_Release_Notes.md#data-release-270) | October 29, 2020 |
| [v26.0](Data_Release_Notes.md#data-release-260) | September 8, 2020 |
| [v25.0](Data_Release_Notes.md#data-release-250) | July 22, 2020 |
| [v24.0](Data_Release_Notes.md#data-release-240) | May 7, 2020 |
| [v23.0](Data_Release_Notes.md#data-release-230) | April 7, 2020 |
| [v22.0](Data_Release_Notes.md#data-release-220) | January 16, 2020 |
| [v21.0](Data_Release_Notes.md#data-release-210) | December 10, 2019 |
| [v20.0](Data_Release_Notes.md#data-release-200) | November 11, 2019 |
| [v19.1](Data_Release_Notes.md#data-release-191) | November 6, 2019 |
| [v19.0](Data_Release_Notes.md#data-release-190) | September 17, 2019 |
| [v18.0](Data_Release_Notes.md#data-release-180) | July 8, 2019 |
| [v17.1](Data_Release_Notes.md#data-release-171) | June 12, 2019 |
| [v17.0](Data_Release_Notes.md#data-release-170) | June 5, 2019 |
| [v16.0](Data_Release_Notes.md#data-release-160) | March 26, 2019 |
| [v15.0](Data_Release_Notes.md#data-release-150) | February 20, 2019 |
| [v14.0](Data_Release_Notes.md#data-release-140) | December 18, 2018 |
| [v13.0](Data_Release_Notes.md#data-release-130) | September 27, 2018 |
| [v12.0](Data_Release_Notes.md#data-release-120) | June 13, 2018 |
| [v11.0](Data_Release_Notes.md#data-release-110) | May 21, 2018 |
| [v10.1](Data_Release_Notes.md#data-release-101) | February 15, 2018 |
| [v10.0](Data_Release_Notes.md#data-release-100) | December 21, 2017 |
| [v9.0](Data_Release_Notes.md#data-release-90) | October 24, 2017 |
| [v8.0](Data_Release_Notes.md#data-release-80) | August 22, 2017 |
| [v7.0](Data_Release_Notes.md#data-release-70) | June 29, 2017 |
| [v6.0](Data_Release_Notes.md#data-release-60) | May 9, 2017 |
| [v5.0](Data_Release_Notes.md#data-release-50) | March 16, 2017 |
| [v4.0](Data_Release_Notes.md#data-release-40) | October 31, 2016 |
| [v3.0](Data_Release_Notes.md#data-release-30) | September 16, 2016 |
| [v2.0](Data_Release_Notes.md#data-release-20) | August 9, 2016 |
| [v1.0](Data_Release_Notes.md#initial-data-release-10) | June 6, 2016 |


## Data Release 38.0

* __GDC Product__: Data
* __Release Date__: August XX, 2023

### New Updates

* New Projects
    * MP2PRT-ALL - Molecular Profiling to Predict Response to Treatment for Acute Lymphoblastic Leukemia - phs002005
        *  1,507 cases
        *  WGS
    * CGCI-HTMCP-DLBCL - HIV+ Tumor Molecular Characterization Project - Diffuse Large B-Cell Lymphoma - phs000235
        *  70 cases
        *  WGS, RNA-Seq, miRNA-Seq, Tissue Slide Images
    * MATCH-B - Genomic Characterization CS-MATCH-0007 Arm B - phs002028
        *  33 cases
        *  WXS, RNA-Seq
    * MATCH-N - Genomic Characterization CS-MATCH-0007 Arm N - phs002151
        *  21 cases
        *  WXS, RNA-Seq

* New Cases from Existing Projects
    * CPTAC-3 - GBM and Kidney cohorts 
    * HCMI-CMDC - 31 cases

* New Data Sets
    * 9,368 WGS alignments from the TCGA program
        * 4,676 Cases
        * 9,368 Aliquots
    * All methylation files that were produced with the SeSAMe pipeline was replaced with a new version.
    * TCGA SNP6 data processed with the ASCAT3 and ABSOLUTE pipelines
    * 172 CEL and birdseed files from TCGA SNP6
    * Release of remaining data for CGCI projects CGCI-BGLSP and CGCI-HTMCP-CC

* New Metadata
    * The `sample_type` field has been broken into four orthogonal fields for easier querying. While the `sample_type` field is still present, we recommend using the following fields: `tissue_type`, `tumor_descriptor`, `specimen_type`, `preservation_method`
    * The QC metrics for applicable BAMs are now queryable through the GDC Data Portal and API.
    * The `msi_status` and `msi_score` fields, which were produced using MSISensor2, are now queryable through the GDC Data Portal and API


A complete list of files included in the GDC Data Portal can be found below:

* [gdc_manifest_2023Aug24_data_release_38.0_active.tsv.gz](gdc_manifest_2023Aug24_data_release_38.0_active.tsv.gz)
* [DR38 Project Level Manifests](DR38_project_manifests.tar.gz)

### Bugs Fixed Since Last Release

* The files produced with the SeSAMe pipeline had unfiltered methylation beta values that should be set as N/A for quality reasons.  These files were replaced.

### Known Issues and Workarounds

* Pathology reports do not have any associated case/biospecimen information in the portal. This information can be found in the reports themselves. <!--SV-2118-->  
* 397 alignments from the TCGA program were found to have contamination values over 0.04 ([alignment list](Contaminated_Alignments.dr32.tsv)). The ensemble MAFs produced by these alignments were removed from the Data Portal.
* One methylation aliquot from the TCGA-COAD project, TCGA-D5-6930-01A-11D-1926-05, was not added to the portal and will be added in a future release.
* The slide image viewer does not display properly for 14 slides, which are identified [here](missing_tiling.txt).  The full slide image can be downloaded as an SVS file.
* The Copy Number Estimate files in GENIE are labeled on the portal as TXT while the files are actually in TSV format.  <!--DAT-2728-->
* Some tumor-only annotated VCFs (not raw VCFs) could have a small proportion of variants that appear twice.  Tumor-only annotated VCFs can be identified by searching for workflow "GATK4 MuTect2 Annotation" <!--SV-1425-->
* The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.
* Some miRNA files with QC failed reads were not swapped in DR11.0.  361 aliquots remain to be swapped in a later release. <!--DAT-1589-->
* Mutation frequency may be underestimated when using MAF files for genes that overlap other genes.  This is because MAF files only record one gene per variant.
* Most intronic mutations are removed for MAF generation.  However, validated variants may rescue these in some cases.  Therefore intronic mutations in MAF files are not representative of those called by mutation callers.
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
* TCGA Projects
    * Incorrect information about treatment may be included for patients within TCGA-HNSC and TCGA-LGG.  Please refer to the clinical XML for accurate information on treatment <!--DAT-2264, DAT-2265-->
    * 74 Diagnostic TCGA slides are attached to a portion rather than a sample like the rest of the diagnostic slides. The reflects how these original samples were handled. <!--SV-1111-->
    * Two tissue slide images are unavailable for download from GDC Data Portal <!--DAT-1439-->
    * The raw and annotated VarScan VCF files for aliquot `TCGA-VR-A8ET-01A-11D-A403-09` are not available. These VCFs files will be replaced in a later release.<!--TT-602, DAT-1489-->
    * Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
    * Tumor grade property is not populated <!--SV-585-->
    * Progression_or_recurrence property is not populated <!--SV-584-->
* TARGET projects
    * 11 bam files for TARGET-NBL RNA-Seq are not available in the GDC Data portal <!--DAT-1476-->
    * There are 5051 TARGET files for which `experimental_strategy`, `data_format`, `platform`, and `data_subtype` are blank <!--SV-944-->
    * There are two cases with identical submitter_id `TARGET-10-PARUYU` <!--SV-940-->
    * Some TARGET cases are missing `days_to_last_follow_up` <!--SV-934-->
    * Some TARGET cases are missing `age_at_diagnosis` <!--SV-933-->
    * Some TARGET files are not connected to all related aliquots <!--SV-929-->
    * Samples of TARGET sample_type `Recurrent Blood Derived Cancer - Bone Marrow` are mislabeled as `Recurrent Blood Derived Cancer - Peripheral Blood`.  A workaround is to look at the sample barcode, which is -04 for `Recurrent Blood Derived Cancer - Bone Marrow`. (e.g. `TARGET-20-PAMYAS-04A-03R`) <!--SV-918-->
    * The latest TARGET data is not yet available at the GDC.  For the complete and latest data, please see the [TARGET Data Matrix](https://ocg.cancer.gov/programs/target/data-matrix).  Data that is not present or is not the most up to date includes:
        *  All microarray data and metadata
        *  All sequencing analyzed data and metadata
        *  1180 of 12063 sequencing runs of raw data
    * Demographic information for some TARGET patients is incorrect.  The correct information can be found in the associated clinical supplement file.  Impacted patients are TARGET-50-PAJNUS. <!--SV-710-->
    * No data from TARGET-MDLS is available.

## Data Release 37.0

* __GDC Product__: Data
* __Release Date__: March 29, 2023

### New Updates

* New Projects
    * APOLLO-LUAD - Proteogenomic characterization of lung adenocarcinoma - phs003011
        *  87 cases
        *  WGS, RNA-Seq
    * CGCI-HTMCP-LC - HIV+ Tumor Molecular Characterization Project - Lung Cancer - phs000530
        *  39 cases
        *  WGS, RNA-Seq, miRNA-Seq, Slide Images
    * MATCH-Q - Genomic Characterization CS-MATCH-0007 Arm Q - phs001926
        *  35 cases
        *  WXS, RNA-Seq
    * MATCH-Y - Genomic Characterization CS-MATCH-0007 Arm Y - phs001904
        *  31 cases
        *  WXS, RNA-Seq

* New Data from Existing Projects
    * CPTAC-3 - 139 new cases and two new snRNA-Seq samples
    * HCMI-CMDC - 118 new cases
    * TCGA-THCA - 941 new WGS alignments
    * TARGET-OS and TARGET-ALL-P2 - Masked Somatic Mutation MAFs are now open access and their mutations now appear in the exploration portal.

* Data Migrated from the Legacy Archive to Active Portal
    * Birdseed files that were generated from Affymetrix SNP6 arrays
    * Additional WGS Alignments are now available for TCGA projects
    * Additional samples from RNA-Seq and WXS are now available for TCGA projects


A complete list of files included in the GDC Data Portal can be found below:

* [gdc_manifest_20230329_data_release_37.0_active.tsv.gz](gdc_manifest_20230329_data_release_37.0_active.tsv.gz)
* [DR37 Project Level Manifests](DR37_project_manifests.tar.gz)

### Unavailable Files

* 56 CPTAC-3 snRNA-Seq files are currently unavailable for download. A list of the affected files can be found [here](CPTAC-3_snRNA-Seq_download_affected_files.txt). These files will be restored for download by the next data release. <!--DAT-3433-->


### Bugs Fixed Since Last Release

* Outcome data for the CPTAC program has been updated.
* The `age_at_index` field was incorrectly reported in days in the GENIE program.  These values have been removed as it contained the same information as the `days_to_birth` field. <!--SV-2032-->


### Known Issues and Workarounds

* The current files produced with the SeSAMe pipeline have unfiltered methylation beta values that should be set as N/A for quality reasons.  These files will be replaced in a future release.
* Pathology reports do not have any associated case/biospecimen information in the portal. This information can be found in the reports themselves. <!--SV-2118-->  
* 397 alignments from the TCGA program were found to have contamination values over 0.04 ([alignment list](Contaminated_Alignments.dr32.tsv)). The ensemble MAFs produced by these alignments were removed from the Data Portal.
* One methylation aliquot from the TCGA-COAD project, TCGA-D5-6930-01A-11D-1926-05, was not added to the portal and will be added in a future release.
* The clinical supplement for TARGET-ALL-P1 is not currently available. It will be made available in a future release.
* The slide image viewer does not display properly for 14 slides, which are identified [here](missing_tiling.txt).  The full slide image can be downloaded as an SVS file.
* The Copy Number Estimate files in GENIE are labeled on the portal as TXT while the files are actually in TSV format.  <!--DAT-2728-->
* Some tumor-only annotated VCFs (not raw VCFs) could have a small proportion of variants that appear twice.  Tumor-only annotated VCFs can be identified by searching for workflow "GATK4 MuTect2 Annotation" <!--SV-1425-->
* The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.
* Some miRNA files with QC failed reads were not swapped in DR11.0.  361 aliquots remain to be swapped in a later release. <!--DAT-1589-->
* Mutation frequency may be underestimated when using MAF files for genes that overlap other genes.  This is because MAF files only record one gene per variant.
* Most intronic mutations are removed for MAF generation.  However, validated variants may rescue these in some cases.  Therefore intronic mutations in MAF files are not representative of those called by mutation callers.
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
* TCGA Projects
    * Incorrect information about treatment may be included for patients within TCGA-HNSC and TCGA-LGG.  Please refer to the clinical XML for accurate information on treatment <!--DAT-2264, DAT-2265-->
    * 74 Diagnostic TCGA slides are attached to a portion rather than a sample like the rest of the diagnostic slides. The reflects how these original samples were handled. <!--SV-1111-->
    * Two tissue slide images are unavailable for download from GDC Data Portal <!--DAT-1439-->
    * The raw and annotated VarScan VCF files for aliquot `TCGA-VR-A8ET-01A-11D-A403-09` are not available. These VCFs files will be replaced in a later release.<!--TT-602, DAT-1489-->
    * Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
    * Tumor grade property is not populated <!--SV-585-->
    * Progression_or_recurrence property is not populated <!--SV-584-->
* TARGET projects
    * TARGET CGI BAMs in the Legacy Archive for the following aliquots should not be used because they were not repaired and concatenated into their original composite BAM files by CGHub.
        * TARGET-20-PASJGZ-04A-02D
        * TARGET-30-PAPTLY-01A-01D
        * TARGET-20-PAEIKD-09A-01D
        * TARGET-20-PASMYS-14A-02D
        * TARGET-20-PAMYAS-14A-02D
        * TARGET-10-PAPZST-09A-01D
    * 11 bam files for TARGET-NBL RNA-Seq are not available in the GDC Data portal <!--DAT-1476-->
    * There are 5051 TARGET files for which `experimental_strategy`, `data_format`, `platform`, and `data_subtype` are blank <!--SV-944-->
    * There are two cases with identical submitter_id `TARGET-10-PARUYU` <!--SV-940-->
    * Some TARGET cases are missing `days_to_last_follow_up` <!--SV-934-->
    * Some TARGET cases are missing `age_at_diagnosis` <!--SV-933-->
    * Some TARGET files are not connected to all related aliquots <!--SV-929-->
    * Samples of TARGET sample_type `Recurrent Blood Derived Cancer - Bone Marrow` are mislabeled as `Recurrent Blood Derived Cancer - Peripheral Blood`.  A workaround is to look at the sample barcode, which is -04 for `Recurrent Blood Derived Cancer - Bone Marrow`. (e.g. `TARGET-20-PAMYAS-04A-03R`) <!--SV-918-->
    * The latest TARGET data is not yet available at the GDC.  For the complete and latest data, please see the [TARGET Data Matrix](https://ocg.cancer.gov/programs/target/data-matrix).  Data that is not present or is not the most up to date includes:
        *  All microarray data and metadata
        *  All sequencing analyzed data and metadata
        *  1180 of 12063 sequencing runs of raw data
    * Demographic information for some TARGET patients is incorrect.  The correct information can be found in the associated clinical supplement file.  Impacted patients are TARGET-50-PAJNUS. <!--SV-710-->
    * No data from TARGET-MDLS is available.
* Issues in the Legacy Archive
    * The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
    * SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
    * Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
    * SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
    * TARGET-MDLS cases do not have disease_type or primary_site populated <!--SV-939-->


## Data Release 36.0

* __GDC Product__: Data
* __Release Date__: December 12, 2022

### New Updates

* New Projects
    * MATCH-Z1D - Genomic Characterization CS-MATCH-0007 Arm Z1D - phs001859
        *  36 cases
        *  WXS, RNA-Seq
    * CDDP_EAGLE-1 - CDDP Integrative Analysis of Lung Adenocarcinoma (Phase 2) - phs001239
        *  50 cases
        *  WXS, WGS, RNA-Seq

* New Data from Existing Projects
    * CMI-MPC - new RNA-Seq and WXS data

* Data Migrated from the Legacy Archive to Active Portal
    * WGS Alignments are now available for 25 TCGA Projects
    * Pathology reports from TCGA
    * Affymetrix SNP6 Genotyping Array CEL files
    * A set of WXS and RNA-Seq samples from TCGA and TARGET that failed harmonization at launch have been rerun and are now available in the active portal.
    * TCGA Bisulfite-Seq files can be downloaded using the following manifests:
        * [TARGET-RT](TARGET-RT.BisulfiteSeq_GDC-Manifest.txt)
        * [TCGA-BLCA](TCGA-BLCA.BisulfiteSeq_GDC-Manifest.txt)
        * [TCGA-BRCA](TCGA-BRCA.BisulfiteSeq_GDC-Manifest.txt)
        * [TCGA-COAD](TCGA-COAD.BisulfiteSeq_GDC-Manifest.txt)
        * [TCGA-GBM](TCGA-GBM.BisulfiteSeq_GDC-Manifest.txt)
        * [TCGA-LUAD](TCGA-LUAD.BisulfiteSeq_GDC-Manifest.txt)
        * [TCGA-LUSC](TCGA-LUSC.BisulfiteSeq_GDC-Manifest.txt)
        * [TCGA-READ](TCGA-READ.BisulfiteSeq_GDC-Manifest.txt)
        * [TCGA-STAD](TCGA-STAD.BisulfiteSeq_GDC-Manifest.txt)
        * [TCGA-UCEC](TCGA-UCEC.BisulfiteSeq_GDC-Manifest.txt)

A complete list of files included in the GDC Data Portal can be found below:

* [gdc_manifest_20221212_data_release_36.0_active.tsv.gz](gdc_manifest_20221212_data_release_36.0_active.tsv.gz)

### Unavailable Files

* None


### Bugs Fixed Since Last Release

* The copy number variation data is now available on the GDC Exploration portal.
* The mutations on GDC Exploration were re-built with the correct gene model.


### Known Issues and Workarounds

* Outcome data for the CPTAC program is not up-to-date. Please visit the [Proteomic Data Commons](https://proteomic.datacommons.cancer.gov/pdc/) for updated outcome data for CPTAC.
* Pathology reports do not have any associated case/biospecimen information in the portal. This information can be found in the reports themselves. <!--SV-2118-->  
* 397 alignments from the TCGA program were found to have contamination values over 0.04 ([alignment list](Contaminated_Alignments.dr32.tsv)). The ensemble MAFs produced by these alignments were removed from the Data Portal.
* One methylation aliquot from the TCGA-COAD project, TCGA-D5-6930-01A-11D-1926-05, was not added to the portal and will be added in a future release.
* The clinical supplement for TARGET-ALL-P1 is not currently available. It will be made available in a future release.
* The slide image viewer does not display properly for 14 slides, which are identified [here](missing_tiling.txt).  The full slide image can be downloaded as an SVS file.
* The Copy Number Estimate files in GENIE are labeled on the portal as TXT while the files are actually in TSV format.  <!--DAT-2728-->
* Some tumor-only annotated VCFs (not raw VCFs) could have a small proportion of variants that appear twice.  Tumor-only annotated VCFs can be identified by searching for workflow "GATK4 MuTect2 Annotation" <!--SV-1425-->
* The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.
* Some miRNA files with QC failed reads were not swapped in DR11.0.  361 aliquots remain to be swapped in a later release. <!--DAT-1589-->
* Mutation frequency may be underestimated when using MAF files for genes that overlap other genes.  This is because MAF files only record one gene per variant.
* Most intronic mutations are removed for MAF generation.  However, validated variants may rescue these in some cases.  Therefore intronic mutations in MAF files are not representative of those called by mutation callers.
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
* TCGA Projects
    * Incorrect information about treatment may be included for patients within TCGA-HNSC and TCGA-LGG.  Please refer to the clinical XML for accurate information on treatment <!--DAT-2264, DAT-2265-->
    * 74 Diagnostic TCGA slides are attached to a portion rather than a sample like the rest of the diagnostic slides. The reflects how these original samples were handled. <!--SV-1111-->
    * Two tissue slide images are unavailable for download from GDC Data Portal <!--DAT-1439-->
    * The raw and annotated VarScan VCF files for aliquot `TCGA-VR-A8ET-01A-11D-A403-09` are not available. These VCFs files will be replaced in a later release.<!--TT-602, DAT-1489-->
    * Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
    * Tumor grade property is not populated <!--SV-585-->
    * Progression_or_recurrence property is not populated <!--SV-584-->
* TARGET projects
    * TARGET CGI BAMs in the Legacy Archive for the following aliquots should not be used because they were not repaired and concatenated into their original composite BAM files by CGHub.
        * TARGET-20-PASJGZ-04A-02D
        * TARGET-30-PAPTLY-01A-01D
        * TARGET-20-PAEIKD-09A-01D
        * TARGET-20-PASMYS-14A-02D
        * TARGET-20-PAMYAS-14A-02D
        * TARGET-10-PAPZST-09A-01D
    * 11 bam files for TARGET-NBL RNA-Seq are not available in the GDC Data portal <!--DAT-1476-->
    * There are 5051 TARGET files for which `experimental_strategy`, `data_format`, `platform`, and `data_subtype` are blank <!--SV-944-->
    * There are two cases with identical submitter_id `TARGET-10-PARUYU` <!--SV-940-->
    * Some TARGET cases are missing `days_to_last_follow_up` <!--SV-934-->
    * Some TARGET cases are missing `age_at_diagnosis` <!--SV-933-->
    * Some TARGET files are not connected to all related aliquots <!--SV-929-->
    * Samples of TARGET sample_type `Recurrent Blood Derived Cancer - Bone Marrow` are mislabeled as `Recurrent Blood Derived Cancer - Peripheral Blood`.  A workaround is to look at the sample barcode, which is -04 for `Recurrent Blood Derived Cancer - Bone Marrow`. (e.g. `TARGET-20-PAMYAS-04A-03R`) <!--SV-918-->
    * The latest TARGET data is not yet available at the GDC.  For the complete and latest data, please see the [TARGET Data Matrix](https://ocg.cancer.gov/programs/target/data-matrix).  Data that is not present or is not the most up to date includes:
        *  All microarray data and metadata
        *  All sequencing analyzed data and metadata
        *  1180 of 12063 sequencing runs of raw data
    * Demographic information for some TARGET patients is incorrect.  The correct information can be found in the associated clinical supplement file.  Impacted patients are TARGET-50-PAJNUS. <!--SV-710-->
    * No data from TARGET-MDLS is available.
* Issues in the Legacy Archive
    * The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
    * SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
    * Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
    * SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
    * TARGET-MDLS cases do not have disease_type or primary_site populated <!--SV-939-->

## Data Release 35.0

* __GDC Product__: Data
* __Release Date__: September 28, 2022

### New Updates

* The SomaticSniper variant calling pipeline was deprecated.  To support this, the following changes were made:
    * All SomaticSniper files no longer appear in the portal, but still can be downloaded using the Data Transfer Tool or API using the original UUID.
    * The aggregated somatic mutation and masked somatic mutation files (multi-caller MAFs) have been replaced to reflect the absence of variants from the SomaticSniper pipeline.
    * The mutations on the exploration portal reflect the above-mentioned masked somatic mutation files.
* 10 snRNA-Seq samples were released from the CPTAC-3 project.
* Additional RNA-Seq samples from 2,082 additional cases are now available for the TARGET-AML project.
* Demographic data has been added for 94 cases in TARGET-ALL-P2 and TARGET-ALL-P3 projects. A list of the updated cases can be found
[here](TARGET_demo_update.dr35.tsv).

A complete list of files included in the GDC Data Portal can be found below:

* [gdc_manifest_20220928_data_release_35.0_active.tsv.gz](gdc_manifest_20220928_data_release_35.0_active.tsv.gz)

### Unavailable Files

* None


### Bugs Fixed Since Last Release

* Data from two HCMI-CMDC aliquots (HCM-BROD-0100-C15-85A-01D-A786-36 and HCM-BROD-0679-C43-85M-01D-A80U-36) were incorrectly selected for inclusion into the Exploration Page in Data Release 32 and has been replaced with the correct aliquots (HCM-BROD-0100-C15-01A-11D-A786-36 and HCM-BROD-0679-C43-06A-11D-A80U-36). <!--DAT-3220-->

### Known Issues and Workarounds

* The mutations on GDC Exploration were built with an incorrect gene model.
    * The mutations are still correct in terms of the gene affected, coordinates, DNA changes, amino acid changes, and impact.
    * Mutations associated with genes that were present in GENCODE v36 and not GENCODE v22 are not displayed. This affects less than 1% of mutations.
    * Files downloaded from the the GDC Repository are not affected by this issue.  This only affects mutations that are downloaded from GDC Exploration.
* Pathology reports do not have any associated case/biospecimen information in the portal. This information can be found in the reports themselves. <!--SV-2118-->  
* 397 alignments from the TCGA program were found to have contamination values over 0.04 ([alignment list](Contaminated_Alignments.dr32.tsv)). The ensemble MAFs produced by these alignments were removed from the Data Portal.
* One methylation aliquot from the TCGA-COAD project, TCGA-D5-6930-01A-11D-1926-05, was not added to the portal and will be added in a future release.
* The clinical supplement for TARGET-ALL-P1 is not currently available. It will be made available in a future release.
* Copy number variations currently do not appear in the Exploration page.  This will be restored in a future release.
* Mutations from SomaticSniper were erroneously labelled as LOH (loss of heterozygosity). This affects the VCF files, MAF files, and may cause SomaticSniper mutations to be absent from ensemble MAFs.
* The slide image viewer does not display properly for 14 slides, which are identified [here](missing_tiling.txt).  The full slide image can be downloaded as an SVS file.
* The Copy Number Estimate files in GENIE are labeled on the portal as TXT while the files are actually in TSV format.  <!--DAT-2728-->
* Some tumor-only annotated VCFs (not raw VCFs) could have a small proportion of variants that appear twice.  Tumor-only annotated VCFs can be identified by searching for workflow "GATK4 MuTect2 Annotation" <!--SV-1425-->
* The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.
* Some miRNA files with QC failed reads were not swapped in DR11.0.  361 aliquots remain to be swapped in a later release. <!--DAT-1589-->
* Mutation frequency may be underestimated when using MAF files for genes that overlap other genes.  This is because MAF files only record one gene per variant.
* Most intronic mutations are removed for MAF generation.  However, validated variants may rescue these in some cases.  Therefore intronic mutations in MAF files are not representative of those called by mutation callers.
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
* TCGA Projects
    * Incorrect information about treatment may be included for patients within TCGA-HNSC and TCGA-LGG.  Please refer to the clinical XML for accurate information on treatment <!--DAT-2264, DAT-2265-->
    * 74 Diagnostic TCGA slides are attached to a portion rather than a sample like the rest of the diagnostic slides. The reflects how these original samples were handled. <!--SV-1111-->
    * Two tissue slide images are unavailable for download from GDC Data Portal <!--DAT-1439-->
    * The raw and annotated VarScan VCF files for aliquot `TCGA-VR-A8ET-01A-11D-A403-09` are not available. These VCFs files will be replaced in a later release.<!--TT-602, DAT-1489-->
    * Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
    * Tumor grade property is not populated <!--SV-585-->
    * Progression_or_recurrence property is not populated <!--SV-584-->
* TARGET projects
    * TARGET CGI BAMs in the Legacy Archive for the following aliquots should not be used because they were not repaired and concatenated into their original composite BAM files by CGHub.
        * TARGET-20-PASJGZ-04A-02D
        * TARGET-30-PAPTLY-01A-01D
        * TARGET-20-PAEIKD-09A-01D
        * TARGET-20-PASMYS-14A-02D
        * TARGET-20-PAMYAS-14A-02D
        * TARGET-10-PAPZST-09A-01D
    * 11 bam files for TARGET-NBL RNA-Seq are not available in the GDC Data portal <!--DAT-1476-->
    * There are 5051 TARGET files for which `experimental_strategy`, `data_format`, `platform`, and `data_subtype` are blank <!--SV-944-->
    * There are two cases with identical submitter_id `TARGET-10-PARUYU` <!--SV-940-->
    * Some TARGET cases are missing `days_to_last_follow_up` <!--SV-934-->
    * Some TARGET cases are missing `age_at_diagnosis` <!--SV-933-->
    * Some TARGET files are not connected to all related aliquots <!--SV-929-->
    * Samples of TARGET sample_type `Recurrent Blood Derived Cancer - Bone Marrow` are mislabeled as `Recurrent Blood Derived Cancer - Peripheral Blood`.  A workaround is to look at the sample barcode, which is -04 for `Recurrent Blood Derived Cancer - Bone Marrow`. (e.g. `TARGET-20-PAMYAS-04A-03R`) <!--SV-918-->
    * The latest TARGET data is not yet available at the GDC.  For the complete and latest data, please see the [TARGET Data Matrix](https://ocg.cancer.gov/programs/target/data-matrix).  Data that is not present or is not the most up to date includes:
        *  All microarray data and metadata
        *  All sequencing analyzed data and metadata
        *  1180 of 12063 sequencing runs of raw data
    * Demographic information for some TARGET patients is incorrect.  The correct information can be found in the associated clinical supplement file.  Impacted patients are TARGET-50-PAJNUS. <!--SV-710-->
    * No data from TARGET-MDLS is available.
* Issues in the Legacy Archive
    * The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
    * SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
    * Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
    * SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
    * TARGET-MDLS cases do not have disease_type or primary_site populated <!--SV-939-->

## Data Release 34.0

* __GDC Product__: Data
* __Release Date__: July 27, 2022

### New updates

* 251 cases from the CPTAC-3 project were added to the portal.  This includes all files associated with these cases.
* 243 cases from the BEATAML1.0-COHORT project were added to the portal.  This includes most of the files associated with these cases.
    * The raw tumor-only VCFs from BEATAML1.0-COHORT are downloadable from the BEATAML1.0-COHORT (2022) publication page [here](https://gdc.cancer.gov/about-data/publications/BEATAML1-0-COHORT-2022) and will be added to the Data Portal in a future release.
* WXS mutations from the BEATAML1.0-COHORT project are now available in the Exploration portal.
* Transcript fusion files are now available for the following projects:
    * BEATAML1.0-COHORT
    * CMI-ASC
    * CMI-MBC
    * CPTAC-2
    * CTSP-DLBCL1
    * MMRF-COMMPASS
    * NCICCR-DLBCL
    * OHSU-CNL
    * ORGANOID-PANCREATIC
    * WCDT-MCRPC

A complete list of files included in the GDC Data Portal can be found below:

* [gdc_manifest_20220727_data_release_34.0_active.tsv.gz](gdc_manifest_20220727_data_release_34.0_active.tsv.gz)

### Unavailable Files

* None


### Bugs Fixed Since Last Release

* Data from two HCMI-CMDC aliquots (HCM-BROD-0100-C15-85A-01D-A786-36 and HCM-BROD-0679-C43-85M-01D-A80U-36) were incorrectly selected for inclusion into the Exploration Page in Data Release 32 and has been replaced with the correct aliquots (HCM-BROD-0100-C15-01A-11D-A786-36 and HCM-BROD-0679-C43-06A-11D-A80U-36). <!--DAT-3220-->

### Known Issues and Workarounds

* Pathology reports do not have any associated case/biospecimen information in the portal. This information can be found in the reports themselves. <!--SV-2118-->  
* 397 alignments from the TCGA program were found to have contamination values over 0.04 ([alignment list](Contaminated_Alignments.dr32.tsv)). The ensemble MAFs produced by these alignments were removed from the Data Portal.
* One methylation aliquot from the TCGA-COAD project, TCGA-D5-6930-01A-11D-1926-05, was not added to the portal and will be added in a future release.
* The clinical supplement for TARGET-ALL-P1 is not currently available. It will be made available in a future release.
* Copy number variations currently do not appear in the Exploration page.  This will be restored in a future release.
* Mutations from SomaticSniper were erroneously labelled as LOH (loss of heterozygosity). This affects the VCF files, MAF files, and may cause SomaticSniper mutations to be absent from ensemble MAFs.
* The slide image viewer does not display properly for 14 slides, which are identified [here](missing_tiling.txt).  The full slide image can be downloaded as an SVS file.
* The Copy Number Estimate files in GENIE are labeled on the portal as TXT while the files are actually in TSV format.  <!--DAT-2728-->
* Some tumor-only annotated VCFs (not raw VCFs) could have a small proportion of variants that appear twice.  Tumor-only annotated VCFs can be identified by searching for workflow "GATK4 MuTect2 Annotation" <!--SV-1425-->
* The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.
* Some miRNA files with QC failed reads were not swapped in DR11.0.  361 aliquots remain to be swapped in a later release. <!--DAT-1589-->
* Mutation frequency may be underestimated when using MAF files for genes that overlap other genes.  This is because MAF files only record one gene per variant.
* Most intronic mutations are removed for MAF generation.  However, validated variants may rescue these in some cases.  Therefore intronic mutations in MAF files are not representative of those called by mutation callers.
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
* TCGA Projects
    * Incorrect information about treatment may be included for patients within TCGA-HNSC and TCGA-LGG.  Please refer to the clinical XML for accurate information on treatment <!--DAT-2264, DAT-2265-->
    * 74 Diagnostic TCGA slides are attached to a portion rather than a sample like the rest of the diagnostic slides. The reflects how these original samples were handled. <!--SV-1111-->
    * Two tissue slide images are unavailable for download from GDC Data Portal <!--DAT-1439-->
    * The raw and annotated VarScan VCF files for aliquot `TCGA-VR-A8ET-01A-11D-A403-09` are not available. These VCFs files will be replaced in a later release.<!--TT-602, DAT-1489-->
    * Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
    * Tumor grade property is not populated <!--SV-585-->
    * Progression_or_recurrence property is not populated <!--SV-584-->
* TARGET projects
    * TARGET CGI BAMs in the Legacy Archive for the following aliquots should not be used because they were not repaired and concatenated into their original composite BAM files by CGHub.
        * TARGET-20-PASJGZ-04A-02D
        * TARGET-30-PAPTLY-01A-01D
        * TARGET-20-PAEIKD-09A-01D
        * TARGET-20-PASMYS-14A-02D
        * TARGET-20-PAMYAS-14A-02D
        * TARGET-10-PAPZST-09A-01D
    * 11 bam files for TARGET-NBL RNA-Seq are not available in the GDC Data portal <!--DAT-1476-->
    * There are 5051 TARGET files for which `experimental_strategy`, `data_format`, `platform`, and `data_subtype` are blank <!--SV-944-->
    * There are two cases with identical submitter_id `TARGET-10-PARUYU` <!--SV-940-->
    * Some TARGET cases are missing `days_to_last_follow_up` <!--SV-934-->
    * Some TARGET cases are missing `age_at_diagnosis` <!--SV-933-->
    * Some TARGET files are not connected to all related aliquots <!--SV-929-->
    * Samples of TARGET sample_type `Recurrent Blood Derived Cancer - Bone Marrow` are mislabeled as `Recurrent Blood Derived Cancer - Peripheral Blood`.  A workaround is to look at the sample barcode, which is -04 for `Recurrent Blood Derived Cancer - Bone Marrow`. (e.g. `TARGET-20-PAMYAS-04A-03R`) <!--SV-918-->
    * The latest TARGET data is not yet available at the GDC.  For the complete and latest data, please see the [TARGET Data Matrix](https://ocg.cancer.gov/programs/target/data-matrix).  Data that is not present or is not the most up to date includes:
        *  All microarray data and metadata
        *  All sequencing analyzed data and metadata
        *  1180 of 12063 sequencing runs of raw data
    * Demographic information for some TARGET patients is incorrect.  The correct information can be found in the associated clinical supplement file.  Impacted patients are TARGET-50-PAJNUS. <!--SV-710-->
    * No data from TARGET-MDLS is available.
* Issues in the Legacy Archive
    * The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
    * SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
    * Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
    * SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
    * TARGET-MDLS cases do not have disease_type or primary_site populated <!--SV-939-->


## Data Release 33.1

* __GDC Product__: Data
* __Release Date__: May 31, 2022

### New updates

* None, see "Bugs Fixed Since Last Release" section below.

A complete list of files included in the GDC Data Portal can be found below:

* [gdc_manifest_20220531_data_release_33.1_active.tsv.gz](gdc_manifest_20220531_data_release_33.1_active.tsv.gz)

### Unavailable Files

* None


### Bugs Fixed Since Last Release

* 32 cases from the EXCEPTIONAL_RESPONDERS-ER project were released as they were missing from the previous release.
* All mutations from EXCEPTIONAL_RESPONDERS-ER in the exploration portal come from WXS data, whereas they were previously a mixture of WXS and Targeted Sequencing.

### Known Issues and Workarounds

* Pathology reports do not have any associated case/biospecimen information in the portal. This information can be found in the reports themselves. <!--SV-2118-->  
* 397 alignments from the TCGA program were found to have contamination values over 0.04 ([alignment list](Contaminated_Alignments.dr32.tsv)). The ensemble MAFs produced by these alignments were removed from the Data Portal.
* One methylation aliquot from the TCGA-COAD project, TCGA-D5-6930-01A-11D-1926-05, was not added to the portal and will be added in a future release.
* The clinical supplement for TARGET-ALL-P1 is not currently available. It will be made available in a future release.
* Copy number variations currently do not appear in the Exploration page.  This will be restored in a future release.
* Mutations from SomaticSniper were erroneously labelled as LOH (loss of heterozygosity). This affects the VCF files, MAF files, and may cause SomaticSniper mutations to be absent from ensemble MAFs.
* The slide image viewer does not display properly for 14 slides, which are identified [here](missing_tiling.txt).  The full slide image can be downloaded as an SVS file.
* The Copy Number Estimate files in GENIE are labeled on the portal as TXT while the files are actually in TSV format.  <!--DAT-2728-->
* Some tumor-only annotated VCFs (not raw VCFs) could have a small proportion of variants that appear twice.  Tumor-only annotated VCFs can be identified by searching for workflow "GATK4 MuTect2 Annotation" <!--SV-1425-->
* The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.
* Some miRNA files with QC failed reads were not swapped in DR11.0.  361 aliquots remain to be swapped in a later release. <!--DAT-1589-->
* Mutation frequency may be underestimated when using MAF files for genes that overlap other genes.  This is because MAF files only record one gene per variant.
* Most intronic mutations are removed for MAF generation.  However, validated variants may rescue these in some cases.  Therefore intronic mutations in MAF files are not representative of those called by mutation callers.
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
* TCGA Projects
    * Incorrect information about treatment may be included for patients within TCGA-HNSC and TCGA-LGG.  Please refer to the clinical XML for accurate information on treatment <!--DAT-2264, DAT-2265-->
    * 74 Diagnostic TCGA slides are attached to a portion rather than a sample like the rest of the diagnostic slides. The reflects how these original samples were handled. <!--SV-1111-->
    * Two tissue slide images are unavailable for download from GDC Data Portal <!--DAT-1439-->
    * The raw and annotated VarScan VCF files for aliquot `TCGA-VR-A8ET-01A-11D-A403-09` are not available. These VCFs files will be replaced in a later release.<!--TT-602, DAT-1489-->
    * Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
    * Tumor grade property is not populated <!--SV-585-->
    * Progression_or_recurrence property is not populated <!--SV-584-->
* TARGET projects
    * TARGET CGI BAMs in the Legacy Archive for the following aliquots should not be used because they were not repaired and concatenated into their original composite BAM files by CGHub.
        * TARGET-20-PASJGZ-04A-02D
        * TARGET-30-PAPTLY-01A-01D
        * TARGET-20-PAEIKD-09A-01D
        * TARGET-20-PASMYS-14A-02D
        * TARGET-20-PAMYAS-14A-02D
        * TARGET-10-PAPZST-09A-01D
    * 11 bam files for TARGET-NBL RNA-Seq are not available in the GDC Data portal <!--DAT-1476-->
    * There are 5051 TARGET files for which `experimental_strategy`, `data_format`, `platform`, and `data_subtype` are blank <!--SV-944-->
    * There are two cases with identical submitter_id `TARGET-10-PARUYU` <!--SV-940-->
    * Some TARGET cases are missing `days_to_last_follow_up` <!--SV-934-->
    * Some TARGET cases are missing `age_at_diagnosis` <!--SV-933-->
    * Some TARGET files are not connected to all related aliquots <!--SV-929-->
    * Samples of TARGET sample_type `Recurrent Blood Derived Cancer - Bone Marrow` are mislabeled as `Recurrent Blood Derived Cancer - Peripheral Blood`.  A workaround is to look at the sample barcode, which is -04 for `Recurrent Blood Derived Cancer - Bone Marrow`. (e.g. `TARGET-20-PAMYAS-04A-03R`) <!--SV-918-->
    * The latest TARGET data is not yet available at the GDC.  For the complete and latest data, please see the [TARGET Data Matrix](https://ocg.cancer.gov/programs/target/data-matrix).  Data that is not present or is not the most up to date includes:
        *  All microarray data and metadata
        *  All sequencing analyzed data and metadata
        *  1180 of 12063 sequencing runs of raw data
    * Demographic information for some TARGET patients is incorrect.  The correct information can be found in the associated clinical supplement file.  Impacted patients are TARGET-50-PAJNUS. <!--SV-710-->
    * No data from TARGET-MDLS is available.
* Issues in the Legacy Archive
    * The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
    * SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
    * Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
    * SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
    * TARGET-MDLS cases do not have disease_type or primary_site populated <!--SV-939-->


## Data Release 33.0

* __GDC Product__: Data
* __Release Date__: May 3, 2022

### New updates

1.  New Project: NCI Exceptional Responders Initiative (EXCEPTIONAL_RESPONDERS-ER, phs001145)
    * RNA-Seq - 45 Cases
    * WXS - 50 Cases
    * Targeted Sequencing - 41 Cases
    * Mutations from WXS and Targeted Sequencing are present in the exploration page.
2.  New Project: Molecular Profiling to Predict Response to Treatment - Wilms Tumor (MP2PRT-WT, phs001965)
    * WGS - 52 Cases
    * RNA-Seq - 52 Cases
    * miRNA-Seq - 52 Cases
3.  Methylation files from the SeSAMe pipeline are now available for CGCI-HTMCP-CC and the TARGET projects.

A complete list of files for this release are listed for the GDC Data Portal and the GDC Legacy Archive are found below:

* [gdc_manifest_20220503_data_release_33.0_active.tsv.gz](gdc_manifest_20220503_data_release_33.0_active.tsv.gz)
* [gdc_manifest_20220503_data_release_33.0_legacy.tsv.gz](gdc_manifest_20220503_data_release_33.0_legacy.tsv.gz)

### Unavailable Files

* The Arriba pipeline failed for one aliquot from EXCEPTIONAL-RESPONDERS-ER and is documented [here](Removed_Aliquots.tsv).


### Bugs Fixed Since Last Release

* Gene-level copy number files from TCGA-THCA and TCGA-UCEC were set as controlled-access files.  These have been corrected to be available as open-access files.
* Due to a problem with the columns generated by the pipeline, all scRNA-Seq files have been replaced with a new version.

### Known Issues and Workarounds

* 397 alignments from the TCGA program were found to have contamination values over 0.04 ([alignment list](Contaminated_Alignments.dr32.tsv)). The ensemble MAFs produced by these alignments were removed from the Data Portal.
* One methylation aliquot from the TCGA-COAD project, TCGA-D5-6930-01A-11D-1926-05, was not added to the portal and will be added in a future release.
* The clinical supplement for TARGET-ALL-P1 is not currently available. It will be made available in a future release.
* Copy number variations currently do not appear in the Exploration page.  This will be restored in a future release.
* Mutations from SomaticSniper were erroneously labelled as LOH (loss of heterozygosity). This affects the VCF files, MAF files, and may cause SomaticSniper mutations to be absent from ensemble MAFs.
* The slide image viewer does not display properly for 14 slides, which are identified [here](missing_tiling.txt).  The full slide image can be downloaded as an SVS file.
* The Copy Number Estimate files in GENIE are labeled on the portal as TXT while the files are actually in TSV format.  <!--DAT-2728-->
* Some tumor-only annotated VCFs (not raw VCFs) could have a small proportion of variants that appear twice.  Tumor-only annotated VCFs can be identified by searching for workflow "GATK4 MuTect2 Annotation" <!--SV-1425-->
* The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.
* Some miRNA files with QC failed reads were not swapped in DR11.0.  361 aliquots remain to be swapped in a later release. <!--DAT-1589-->
* Mutation frequency may be underestimated when using MAF files for genes that overlap other genes.  This is because MAF files only record one gene per variant.
* Most intronic mutations are removed for MAF generation.  However, validated variants may rescue these in some cases.  Therefore intronic mutations in MAF files are not representative of those called by mutation callers.
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
* TCGA Projects
    * Incorrect information about treatment may be included for patients within TCGA-HNSC and TCGA-LGG.  Please refer to the clinical XML for accurate information on treatment <!--DAT-2264, DAT-2265-->
    * 74 Diagnostic TCGA slides are attached to a portion rather than a sample like the rest of the diagnostic slides. The reflects how these original samples were handled. <!--SV-1111-->
    * Two tissue slide images are unavailable for download from GDC Data Portal <!--DAT-1439-->
    * The raw and annotated VarScan VCF files for aliquot `TCGA-VR-A8ET-01A-11D-A403-09` are not available. These VCFs files will be replaced in a later release.<!--TT-602, DAT-1489-->
    * Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
    * Tumor grade property is not populated <!--SV-585-->
    * Progression_or_recurrence property is not populated <!--SV-584-->
* TARGET projects
    * TARGET CGI BAMs in the Legacy Archive for the following aliquots should not be used because they were not repaired and concatenated into their original composite BAM files by CGHub.
        * TARGET-20-PASJGZ-04A-02D
        * TARGET-30-PAPTLY-01A-01D
        * TARGET-20-PAEIKD-09A-01D
        * TARGET-20-PASMYS-14A-02D
        * TARGET-20-PAMYAS-14A-02D
        * TARGET-10-PAPZST-09A-01D
    * 11 bam files for TARGET-NBL RNA-Seq are not available in the GDC Data portal <!--DAT-1476-->
    * There are 5051 TARGET files for which `experimental_strategy`, `data_format`, `platform`, and `data_subtype` are blank <!--SV-944-->
    * There are two cases with identical submitter_id `TARGET-10-PARUYU` <!--SV-940-->
    * Some TARGET cases are missing `days_to_last_follow_up` <!--SV-934-->
    * Some TARGET cases are missing `age_at_diagnosis` <!--SV-933-->
    * Some TARGET files are not connected to all related aliquots <!--SV-929-->
    * Samples of TARGET sample_type `Recurrent Blood Derived Cancer - Bone Marrow` are mislabeled as `Recurrent Blood Derived Cancer - Peripheral Blood`.  A workaround is to look at the sample barcode, which is -04 for `Recurrent Blood Derived Cancer - Bone Marrow`. (e.g. `TARGET-20-PAMYAS-04A-03R`) <!--SV-918-->
    * The latest TARGET data is not yet available at the GDC.  For the complete and latest data, please see the [TARGET Data Matrix](https://ocg.cancer.gov/programs/target/data-matrix).  Data that is not present or is not the most up to date includes:
        *  All microarray data and metadata
        *  All sequencing analyzed data and metadata
        *  1180 of 12063 sequencing runs of raw data
    * Demographic information for some TARGET patients is incorrect.  The correct information can be found in the associated clinical supplement file.  Impacted patients are TARGET-50-PAJNUS. <!--SV-710-->
    * No data from TARGET-MDLS is available.
* Issues in the Legacy Archive
    * The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
    * SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
    * Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
    * SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
    * TARGET-MDLS cases do not have disease_type or primary_site populated <!--SV-939-->

## Data Release 32.0

* __GDC Product__: Data - GENCODE v36 Release
* __Release Date__: March 29, 2022

### New updates

#### New data files
1.  The following data types have been replaced with new GENCODE v36 versions
    * RNA-Seq: all files, including alignments, gene expression files, and transcript fusion files.
    * WXS and Targeted Sequencing: annotated VCFs, single-caller MAFs, Ensemble MAFs.
    * WGS: BEDPE-format structural variants and gene-level copy number variants.
    * GENIE Targeted Sequencing files.
    * FM-AD Targeted Sequencing files.
        * The primary-site-level FM-AD MAF files have been replaced with aliquot-level MAF files.
1. RNA-Seq STAR-Counts files now contain additional normalized counts such as FPKM, FPKM-UQ, and TPM.
1. All WXS files for TCGA have been replaced with new versions. Alignments will contain QC metrics and variants were produced using the same pipelines as all other GDC projects.
1. TCGA RNA-Seq has been changed to contain three alignments (genomic, transcriptome, and chimeric), STAR-counts files, and transcript fusion files for each aliquot.
1. The project-level MAFs in TCGA and FM-AD have been replaced with aliquot-level MAFs.
1. GENCODE v22 derived files (not BAM) that no longer appear in the portal will be downloadable as previous versions of v36 files.  
1. Methylation data produced from the SeSAMe pipeline is now available for all TCGA projects.
1. Note that miRNA-Seq data remains unchanged. The miRNA-Seq pipeline uses the miRBase database, which is not affected by the GENCODE version change.
1. A set of manifests were generated at the project-level that map each v22 file to its corresponding v36 file.  These can be used to help users transition from v22 to v36 and can be downloaded [here](https://github.com/NCI-GDC/gdc-docs/tree/develop/docs/Data/Release_Notes/GCv36_Manifests).

#### Removed data files and pipelines
1. Files from the HTSeq pipeline are no longer supported and will no longer appear in the portal. Normalized counts can now be found in the STAR-Counts files.  
1. Files that originated from the methylation liftover pipeline are no longer supported and will no longer appear in the portal.
1. GENCODE v22 BAM files that no longer appear in the portal will be available for six months past this release. They may not be available after that.
1. New variant calling tumor-normal pairing was implemented in TCGA, which results in certain aliquots no longer being available as a v36 version (see the aliquots labeled "Unpaired Aliquots" [here](Removed_Aliquots.tsv)).
1. Some aliquots failed harmonization when the new v36 gene model was used, which results in some new versions no longer being available (see the aliquots labeled "Failed Harmonization" [here](Removed_Aliquots.tsv)).
1. Some aliquots were found to contain a cross-patient contamination level of over 0.04 as measured by GATK4 CalculateContamination (see the aliquots labeled "Contamination" [here](Removed_Aliquots.tsv)).  

#### Data Portal Exploration Data
1. The Data Portal Exploration Page is now populated based on open-access mutations from analyses that used GENCODE v36.
1. Mutations from SomaticSniper will not appear on the Exploration page.
1. Due to the copy number variation pipeline transition from GISTIC to ASCAT, the CNV data was not included in the GDC Exploration page. This will be replaced in a future release once visualization of the new pipeline is fully assessed.
1. The TCGA program mutations have been processed using the same pipeline as all other projects, which resulted in a 26% reduction in the number of open-access mutations. Some points on this change are listed below with TCGA-BRCA as the benchmark project:
    * 97% of the previously released open-access mutations are still discoverable in the new GDC controlled-access MAFs. This number increases to 99.95% when focusing only on mutations that were also called by [MC3](https://gdc.cancer.gov/about-data/publications/pancanatlas).
    * Somatic mutations will now be removed from the Data Portal Exploration Page unless they are detected by more than one variant calling software. This accounts for 40% of the total reduction.
    * Somatic mutations will now be removed from the Data Portal Exploration Page if they are detected outside of the target capture region, while previously out-of-target mutations detected from the TCGA Gene Annotation File (GAF) regions were allowed. This accounts for 36% of the total reduction.  
    * Some TCGA-specific variant-rescue steps have been removed in favor of a more robust and uniform filtering pipeline.
    * Some other minor changes due to updates in the gene model or other databases (e.g., the ExAC germline variant database was replaced with gnomAD in DR32).


A complete list of files for this release are listed for the GDC Data Portal and the GDC Legacy Archive are found below:

* [gdc_manifest_20220316_data_release_32.0_active.tsv.gz](gdc_manifest_20220316_data_release_32.0_active.tsv.gz)
* [gdc_manifest_20220316_data_release_32.0_legacy.tsv.gz](gdc_manifest_20220316_data_release_32.0_legacy.tsv.gz)

### Bugs Fixed Since Last Release

* None

### Known Issues and Workarounds

* 397 alignments from the TCGA program were found to have contamination values over 0.04 ([alignment list](Contaminated_Alignments.dr32.tsv)). The ensemble MAFs produced by these alignments were removed from the Data Portal.
* One methylation aliquot from the TCGA-COAD project, TCGA-D5-6930-01A-11D-1926-05, was not added to the portal and will be added in a future release.
* The clinical supplement for TARGET-ALL-P1 is not currently available. It will be made available in a future release.
* Copy number variations currently do not appear in the Exploration page.  This will be restored in a future release.
* Mutations from SomaticSniper were erroneously labelled as LOH (loss of heterozygosity). This affects the VCF files, MAF files, and may cause SomaticSniper mutations to be absent from ensemble MAFs.
* The slide image viewer does not display properly for 14 slides, which are identified [here](missing_tiling.txt).  The full slide image can be downloaded as an SVS file.
* The Copy Number Estimate files in GENIE are labeled on the portal as TXT while the files are actually in TSV format.  <!--DAT-2728-->
* Some tumor-only annotated VCFs (not raw VCFs) could have a small proportion of variants that appear twice.  Tumor-only annotated VCFs can be identified by searching for workflow "GATK4 MuTect2 Annotation" <!--SV-1425-->
* The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.
* Some miRNA files with QC failed reads were not swapped in DR11.0.  361 aliquots remain to be swapped in a later release. <!--DAT-1589-->
* Mutation frequency may be underestimated when using MAF files for genes that overlap other genes.  This is because MAF files only record one gene per variant.
* Most intronic mutations are removed for MAF generation.  However, validated variants may rescue these in some cases.  Therefore intronic mutations in MAF files are not representative of those called by mutation callers.
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
* TCGA Projects
    * Incorrect information about treatment may be included for patients within TCGA-HNSC and TCGA-LGG.  Please refer to the clinical XML for accurate information on treatment <!--DAT-2264, DAT-2265-->
    * 74 Diagnostic TCGA slides are attached to a portion rather than a sample like the rest of the diagnostic slides. The reflects how these original samples were handled. <!--SV-1111-->
    * Two tissue slide images are unavailable for download from GDC Data Portal <!--DAT-1439-->
    * The raw and annotated VarScan VCF files for aliquot `TCGA-VR-A8ET-01A-11D-A403-09` are not available. These VCFs files will be replaced in a later release.<!--TT-602, DAT-1489-->
    * Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
    * Tumor grade property is not populated <!--SV-585-->
    * Progression_or_recurrence property is not populated <!--SV-584-->
* TARGET projects
    * TARGET CGI BAMs in the Legacy Archive for the following aliquots should not be used because they were not repaired and concatenated into their original composite BAM files by CGHub.
        * TARGET-20-PASJGZ-04A-02D
        * TARGET-30-PAPTLY-01A-01D
        * TARGET-20-PAEIKD-09A-01D
        * TARGET-20-PASMYS-14A-02D
        * TARGET-20-PAMYAS-14A-02D
        * TARGET-10-PAPZST-09A-01D
    * 11 bam files for TARGET-NBL RNA-Seq are not available in the GDC Data portal <!--DAT-1476-->
    * There are 5051 TARGET files for which `experimental_strategy`, `data_format`, `platform`, and `data_subtype` are blank <!--SV-944-->
    * There are two cases with identical submitter_id `TARGET-10-PARUYU` <!--SV-940-->
    * Some TARGET cases are missing `days_to_last_follow_up` <!--SV-934-->
    * Some TARGET cases are missing `age_at_diagnosis` <!--SV-933-->
    * Some TARGET files are not connected to all related aliquots <!--SV-929-->
    * Samples of TARGET sample_type `Recurrent Blood Derived Cancer - Bone Marrow` are mislabeled as `Recurrent Blood Derived Cancer - Peripheral Blood`.  A workaround is to look at the sample barcode, which is -04 for `Recurrent Blood Derived Cancer - Bone Marrow`. (e.g. `TARGET-20-PAMYAS-04A-03R`) <!--SV-918-->
    * The latest TARGET data is not yet available at the GDC.  For the complete and latest data, please see the [TARGET Data Matrix](https://ocg.cancer.gov/programs/target/data-matrix).  Data that is not present or is not the most up to date includes:
        *  All microarray data and metadata
        *  All sequencing analyzed data and metadata
        *  1180 of 12063 sequencing runs of raw data
    * Demographic information for some TARGET patients is incorrect.  The correct information can be found in the associated clinical supplement file.  Impacted patients are TARGET-50-PAJNUS. <!--SV-710-->
    * No data from TARGET-MDLS is available.
* Issues in the Legacy Archive
    * The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
    * SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
    * Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
    * SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
    * TARGET-MDLS cases do not have disease_type or primary_site populated <!--SV-939-->




## Data Release 31.0

* __GDC Product__: Data
* __Release Date__: October 29, 2021

### New updates

1.  TCGA Slide Images:
    * All TCGA slide images that were removed earlier this year have been restored.
    * Note that the UUIDs for most TCGA slide images have changed.  Older manifest files may not work when downloading slide images.
2. CPTAC-3 clinical data has been refreshed and includes new follow up entities.
3. REBC-THYR
    * The clinical and biospecimen XML files were removed as they were not intended for release in DR 30.
    * The case REBC-ADL5 was added, which includes one WGS pair.

A complete list of files for this release are listed for the GDC Data Portal and the GDC Legacy Archive are found below:

* [gdc_manifest_20211029_data_release_31.0_active.tsv.gz](gdc_manifest_20211029_data_release_31.0_active.tsv.gz)
* [gdc_manifest_20211029_data_release_31.0_legacy.tsv.gz](gdc_manifest_20211029_data_release_31.0_legacy.tsv.gz)

### Bugs Fixed Since Last Release

* One file from a previous version of the methylation pipeline appeared in the data portal (bd2f864a-3f00-47b5-815d-bd01ca21ef61; CPTAC-3).  This file should no longer appear in the data portal. <!--SV-2012-->

### Known Issues and Workarounds

* The slide image viewer does not display properly for 14 slides, which are identified [here](missing_tiling.txt).  The full slide image can be downloaded as an SVS file.
* The Copy Number Estimate files in GENIE are labeled on the portal as TXT while the files are actually in TSV format.  <!--DAT-2728-->
* Some tumor-only annotated VCFs (not raw VCFs) could have a small proportion of variants that appear twice.  Tumor-only annotated VCFs can be identified by searching for workflow "GATK4 MuTect2 Annotation" <!--SV-1425-->
* The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.
* Some miRNA files with QC failed reads were not swapped in DR11.0.  361 aliquots remain to be swapped in a later release <!--DAT-1589-->
* Mutation frequency may be underestimated when using MAF files for genes that overlap other genes.  This is because MAF files only record one gene per variant.
* Most intronic mutations are removed for MAF generation.  However, validated variants may rescue these in some cases.  Therefore intronic mutations in MAF files are not representative of those called by mutation callers.
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
* TCGA Projects
    * Incorrect information about treatment may be included for patients within TCGA-HNSC and TCGA-LGG.  Please refer to the clinical XML for accurate information on treatment <!--DAT-2264, DAT-2265-->
    * 74 Diagnostic TCGA slides are attached to a portion rather than a sample like the rest of the diagnostic slides. The reflects how these original samples were handled. <!--SV-1111-->
    * Two tissue slide images are unavailable for download from GDC Data Portal <!--DAT-1439-->
    * The raw and annotated VarScan VCF files for aliquot `TCGA-VR-A8ET-01A-11D-A403-09` are not available. These VCFs files will be replaced in a later release.<!--TT-602, DAT-1489-->
    * Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
    * Tumor grade property is not populated <!--SV-585-->
    * Progression_or_recurrence property is not populated <!--SV-584-->
* TARGET projects
    * TARGET CGI BAMs in the Legacy Archive for the following aliquots should not be used because they were not repaired and concatenated into their original composite BAM files by CGHub.
        * TARGET-20-PASJGZ-04A-02D
        * TARGET-30-PAPTLY-01A-01D
        * TARGET-20-PAEIKD-09A-01D
        * TARGET-20-PASMYS-14A-02D
        * TARGET-20-PAMYAS-14A-02D
        * TARGET-10-PAPZST-09A-01D
    * 11 bam files for TARGET-NBL RNA-Seq are not available in the GDC Data portal <!--DAT-1476-->
    * There are 5051 TARGET files for which `experimental_strategy`, `data_format`, `platform`, and `data_subtype` are blank <!--SV-944-->
    * There are two cases with identical submitter_id `TARGET-10-PARUYU` <!--SV-940-->
    * Some TARGET cases are missing `days_to_last_follow_up` <!--SV-934-->
    * Some TARGET cases are missing `age_at_diagnosis` <!--SV-933-->
    * Some TARGET files are not connected to all related aliquots <!--SV-929-->
    * Samples of TARGET sample_type `Recurrent Blood Derived Cancer - Bone Marrow` are mislabeled as `Recurrent Blood Derived Cancer - Peripheral Blood`.  A workaround is to look at the sample barcode, which is -04 for `Recurrent Blood Derived Cancer - Bone Marrow`. (e.g. `TARGET-20-PAMYAS-04A-03R`) <!--SV-918-->
    * The latest TARGET data is not yet available at the GDC.  For the complete and latest data, please see the [TARGET Data Matrix](https://ocg.cancer.gov/programs/target/data-matrix).  Data that is not present or is not the most up to date includes:
        *  All microarray data and metadata
        *  All sequencing analyzed data and metadata
        *  1180 of 12063 sequencing runs of raw data
    * Demographic information for some TARGET patients is incorrect.  The correct information can be found in the associated clinical supplement file.  Impacted patients are TARGET-50-PAJNUS. <!--SV-710-->
    * No data from TARGET-MDLS is available.
* Issues in the Legacy Archive
    * The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
    * SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
    * Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
    * SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
    * TARGET-MDLS cases do not have disease_type or primary_site populated <!--SV-939-->

## Data Release 30.0

* __GDC Product__: Data
* __Release Date__: September 23, 2021

### New updates

1.  New Projects:
    * TRIO-CRU (phs001163) - Ukrainian National Research Center for Radiation Medicine Trio Study
        * WGS Alignments
    * REBC-THYR (phs001134) - Comprehensive genomic characterization of radiation-related papillary thyroid cancer in the Ukraine
        * miRNA-Seq
        * RNA-Seq
        * WGS
2.  CPTAC Program
    * CPTAC-3 methylation data produced from the SeSAMe pipeline is now available.
    * CPTAC-2 miRNA-Seq files have been replaced with better quality data.
3. HCMI-CMDC
    * 31 New cases have been released to the GDC Data Portal.
    * Methylation data produced from the SeSAMe pipeline is now available.
4. TCGA
    * Protein expression data (RPPA) is now available for 32 projects.
    * RNA-Seq data for TCGA-TGCT was replaced with files from an updated pipeline.
5. TARGET-AML - New RNA-Seq and miRNA-Seq aliquots have been released.

A complete list of files for this release are listed for the GDC Data Portal and the GDC Legacy Archive are found below:

* [gdc_manifest_20210923_data_release_30.0_active.tsv.gz](gdc_manifest_20210923_data_release_30.0_active.tsv.gz)
* [gdc_manifest_20210923_data_release_30.0_legacy.tsv.gz](gdc_manifest_20210923_data_release_30.0_legacy.tsv.gz)

### Bugs Fixed Since Last Release

*  None

### Known Issues and Workarounds

* One file from a previous version of the methylation pipeline appears in the data portal (bd2f864a-3f00-47b5-815d-bd01ca21ef61; CPTAC-3).  This file cannot be downloaded, but may cause bulk downloads to fail.  Remove this file from any manifest or cart you plan on downloading. <!--SV-2012-->
* The Copy Number Estimate files in GENIE are labeled on the portal as TXT while the files are actually in TSV format.  <!--DAT-2728-->
* Some tumor-only annotated VCFs (not raw VCFs) could have a small proportion of variants that appear twice.  Tumor-only annotated VCFs can be identified by searching for workflow "GATK4 MuTect2 Annotation" <!--SV-1425-->
* The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.
* Some miRNA files with QC failed reads were not swapped in DR11.0.  361 aliquots remain to be swapped in a later release <!--DAT-1589-->
* Mutation frequency may be underestimated when using MAF files for genes that overlap other genes.  This is because MAF files only record one gene per variant.
* Most intronic mutations are removed for MAF generation.  However, validated variants may rescue these in some cases.  Therefore intronic mutations in MAF files are not representative of those called by mutation callers.
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
* TCGA Projects
    * Incorrect information about treatment may be included for patients within TCGA-HNSC and TCGA-LGG.  Please refer to the clinical XML for accurate information on treatment <!--DAT-2264, DAT-2265-->
    * 74 Diagnostic TCGA slides are attached to a portion rather than a sample like the rest of the diagnostic slides. The reflects how these original samples were handled. <!--SV-1111-->
    * Two tissue slide images are unavailable for download from GDC Data Portal <!--DAT-1439-->
    * The raw and annotated VarScan VCF files for aliquot `TCGA-VR-A8ET-01A-11D-A403-09` are not available. These VCFs files will be replaced in a later release.<!--TT-602, DAT-1489-->
    * Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
    * Tumor grade property is not populated <!--SV-585-->
    * Progression_or_recurrence property is not populated <!--SV-584-->
* TARGET projects
    * TARGET CGI BAMs in the Legacy Archive for the following aliquots should not be used because they were not repaired and concatenated into their original composite BAM files by CGHub.
        * TARGET-20-PASJGZ-04A-02D
        * TARGET-30-PAPTLY-01A-01D
        * TARGET-20-PAEIKD-09A-01D
        * TARGET-20-PASMYS-14A-02D
        * TARGET-20-PAMYAS-14A-02D
        * TARGET-10-PAPZST-09A-01D
    * 11 bam files for TARGET-NBL RNA-Seq are not available in the GDC Data portal <!--DAT-1476-->
    * There are 5051 TARGET files for which `experimental_strategy`, `data_format`, `platform`, and `data_subtype` are blank <!--SV-944-->
    * There are two cases with identical submitter_id `TARGET-10-PARUYU` <!--SV-940-->
    * Some TARGET cases are missing `days_to_last_follow_up` <!--SV-934-->
    * Some TARGET cases are missing `age_at_diagnosis` <!--SV-933-->
    * Some TARGET files are not connected to all related aliquots <!--SV-929-->
    * Samples of TARGET sample_type `Recurrent Blood Derived Cancer - Bone Marrow` are mislabeled as `Recurrent Blood Derived Cancer - Peripheral Blood`.  A workaround is to look at the sample barcode, which is -04 for `Recurrent Blood Derived Cancer - Bone Marrow`. (e.g. `TARGET-20-PAMYAS-04A-03R`) <!--SV-918-->
    * The latest TARGET data is not yet available at the GDC.  For the complete and latest data, please see the [TARGET Data Matrix](https://ocg.cancer.gov/programs/target/data-matrix).  Data that is not present or is not the most up to date includes:
        *  All microarray data and metadata
        *  All sequencing analyzed data and metadata
        *  1180 of 12063 sequencing runs of raw data
    * Demographic information for some TARGET patients is incorrect.  The correct information can be found in the associated clinical supplement file.  Impacted patients are TARGET-50-PAJNUS. <!--SV-710-->
    * No data from TARGET-MDLS is available.
* Issues in the Legacy Archive
    * The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
    * SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
    * Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
    * SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
    * TARGET-MDLS cases do not have disease_type or primary_site populated <!--SV-939-->


## Data Release 29.0

* __GDC Product__: Data
* __Release Date__: March 31, 2021

### New updates

1.  Count Me In Program
    * Aliquot-level MAFs are now available for projects CMI-ASC, CMI-MBC, and CMI-MPC.
    * Somatic mutation are now explorable for projects CMI-ASC, CMI-MBC, and CMI-MPC
2.  CPTAC Program
    * CPTAC-2 open-access somatic mutations are now browsable through the GDC Exploration Portal.
    * MSI data is now browsable through the faceted search for CPTAC-2 and CPTAC-3.
3. HCMI-CMDC - Data files and explorable mutations for 18 new cases are now available.

A complete list of files for this release are listed for the GDC Data Portal and the GDC Legacy Archive are found below:

* [gdc_manifest_20210331_data_release_29.0_active.tsv.gz](gdc_manifest_20210331_data_release_29.0_active.tsv.gz)
* [gdc_manifest_20210331_data_release_29.0_legacy.tsv.gz](gdc_manifest_20210331_data_release_29.0_legacy.tsv.gz)

### Bugs Fixed Since Last Release

*  The aggregated and masked MAF files that were missing for seven pancreatic cases in CPTAC-3 have been restored to the data portal.
* The missing RNA-Seq data files for the seven normal pancreatic cases in CPTAC-3 have been restored to the data portal.

### Known Issues and Workarounds

* The Copy Number Estimate files in GENIE are labeled on the portal as TXT while the files are actually in TSV format.  <!--DAT-2728-->
* Some tumor-only annotated VCFs (not raw VCFs) could have a small proportion of variants that appear twice.  Tumor-only annotated VCFs can be identified by searching for workflow "GATK4 MuTect2 Annotation" <!--SV-1425-->
* The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.
* Some miRNA files with QC failed reads were not swapped in DR11.0.  361 aliquots remain to be swapped in a later release <!--DAT-1589-->
* Mutation frequency may be underestimated when using MAF files for genes that overlap other genes.  This is because MAF files only record one gene per variant.
* Most intronic mutations are removed for MAF generation.  However, validated variants may rescue these in some cases.  Therefore intronic mutations in MAF files are not representative of those called by mutation callers.
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
* TCGA Projects
    * Incorrect information about treatment may be included for patients within TCGA-HNSC and TCGA-LGG.  Please refer to the clinical XML for accurate information on treatment <!--DAT-2264, DAT-2265-->
    * 74 Diagnostic TCGA slides are attached to a portion rather than a sample like the rest of the diagnostic slides. The reflects how these original samples were handled. <!--SV-1111-->
    * Two tissue slide images are unavailable for download from GDC Data Portal <!--DAT-1439-->
    * The raw and annotated VarScan VCF files for aliquot `TCGA-VR-A8ET-01A-11D-A403-09` are not available. These VCFs files will be replaced in a later release.<!--TT-602, DAT-1489-->
    * Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
    * Tumor grade property is not populated <!--SV-585-->
    * Progression_or_recurrence property is not populated <!--SV-584-->
* TARGET projects
    * TARGET CGI BAMs in the Legacy Archive for the following aliquots should not be used because they were not repaired and concatenated into their original composite BAM files by CGHub.
        * TARGET-20-PASJGZ-04A-02D
        * TARGET-30-PAPTLY-01A-01D
        * TARGET-20-PAEIKD-09A-01D
        * TARGET-20-PASMYS-14A-02D
        * TARGET-20-PAMYAS-14A-02D
        * TARGET-10-PAPZST-09A-01D
    * 11 bam files for TARGET-NBL RNA-Seq are not available in the GDC Data portal <!--DAT-1476-->
    * There are 5051 TARGET files for which `experimental_strategy`, `data_format`, `platform`, and `data_subtype` are blank <!--SV-944-->
    * There are two cases with identical submitter_id `TARGET-10-PARUYU` <!--SV-940-->
    * Some TARGET cases are missing `days_to_last_follow_up` <!--SV-934-->
    * Some TARGET cases are missing `age_at_diagnosis` <!--SV-933-->
    * Some TARGET files are not connected to all related aliquots <!--SV-929-->
    * Samples of TARGET sample_type `Recurrent Blood Derived Cancer - Bone Marrow` are mislabeled as `Recurrent Blood Derived Cancer - Peripheral Blood`.  A workaround is to look at the sample barcode, which is -04 for `Recurrent Blood Derived Cancer - Bone Marrow`. (e.g. `TARGET-20-PAMYAS-04A-03R`) <!--SV-918-->
    * The latest TARGET data is not yet available at the GDC.  For the complete and latest data, please see the [TARGET Data Matrix](https://ocg.cancer.gov/programs/target/data-matrix).  Data that is not present or is not the most up to date includes:
        *  All microarray data and metadata
        *  All sequencing analyzed data and metadata
        *  1180 of 12063 sequencing runs of raw data
    * Demographic information for some TARGET patients is incorrect.  The correct information can be found in the associated clinical supplement file.  Impacted patients are TARGET-50-PAJNUS. <!--SV-710-->
    * No data from TARGET-MDLS is available.
* Issues in the Legacy Archive
    * The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
    * SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
    * Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
    * SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
    * TARGET-MDLS cases do not have disease_type or primary_site populated <!--SV-939-->

    ## Data Release 28.0

    * __GDC Product__: Data
    * __Release Date__: February 2, 2021

    ### New updates

    1.  New Project: CMI-MPC - Count Me In - The Metastatic Prostate Cancer Project
        * WXS alignments and variant calls (VCFs) are available.
    2.  New Data Type: Single nuclei (snRNA-Seq) data is now available for 18 CPTAC-3 cases. See the [RNA-Seq](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/#scrna-seq-pipeline) documentation for details.
    3.  CPTAC-3
        * Data files for 147 new cases from the pancreatic cohort are now available.
        * CPTAC-3 open-access somatic mutations are now browsable through the GDC Exploration Portal.
        * RNA-Seq transcript fusion files are now available.
        * Targeted Sequencing alignments and raw tumor-only variant calls (VCF) are now available.
    4. HCMI-CMDC
        * Data files for 22 new cases are now available.
        * The HCMI-CMDC open-access somatic mutations have been refreshed on the GDC Exploration Portal to reflect all newly released cases.

    A complete list of files for this release are listed for the GDC Data Portal and the GDC Legacy Archive are found below:

    * [gdc_manifest_20210202_data_release_28.0_active.tsv.gz](gdc_manifest_20210202_data_release_28.0_active.tsv.gz)
    * [gdc_manifest_20210202_data_release_28.0_legacy.tsv.gz](gdc_manifest_20210202_data_release_28.0_legacy.tsv.gz)

    ### Bugs Fixed Since Last Release

    *  None

    ### Known Issues and Workarounds

    * The aggregated and masked MAF files for seven pancreatic cases in CPTAC-3 do not appear in the Data Portal. See below for download instructions.
        - [This manifest](CPTAC-3_7Cases-WXS-MAFs_GDC-Manifest.txt) can be used to download the files.  
        - To download the raw aggregated MAF files, dbGaP access to CPTAC-3 (phs001287) is required.  The masked MAF files are open-access.
        - The seven cases are as follows: C3L-04027, C3L-04080, C3N-02585, C3N-02768, C3N-02971, C3N-03754, and C3N-03839. The case the each file is associated with is denoted in the manifest.
    * The RNA-Seq data files for the seven normal pancreatic cases in CPTAC-3 do not appear in the Data Portal. See below for download instructions.
        - [This manifest](CPTAC-3_7CasesRNASeq_GDC-Manifest.txt) can be used to download the files.  
        - To download the alignments or splice-junction files, dbGaP access to CPTAC-3 (phs001287) is required.  The other gene expression files are open-access.
        - The seven cases are as follows: C3L-03513, C3L-07032, C3L-07033, C3L-07034, C3L-07035, C3L-07036, C3L-07037. The case the each file is associated with is denoted in the manifest.
    * The Copy Number Estimate files in GENIE are labeled on the portal as TXT while the files are actually in TSV format.  <!--DAT-2728-->
    * Some tumor-only annotated VCFs (not raw VCFs) could have a small proportion of variants that appear twice.  Tumor-only annotated VCFs can be identified by searching for workflow "GATK4 MuTect2 Annotation" <!--SV-1425-->
    * The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.
    * Some miRNA files with QC failed reads were not swapped in DR11.0.  361 aliquots remain to be swapped in a later release <!--DAT-1589-->
    * Mutation frequency may be underestimated when using MAF files for genes that overlap other genes.  This is because MAF files only record one gene per variant.
    * Most intronic mutations are removed for MAF generation.  However, validated variants may rescue these in some cases.  Therefore intronic mutations in MAF files are not representative of those called by mutation callers.
    * Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
    * BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
    * Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
    * TCGA Projects
        * Incorrect information about treatment may be included for patients within TCGA-HNSC and TCGA-LGG.  Please refer to the clinical XML for accurate information on treatment <!--DAT-2264, DAT-2265-->
        * 74 Diagnostic TCGA slides are attached to a portion rather than a sample like the rest of the diagnostic slides. The reflects how these original samples were handled. <!--SV-1111-->
        * Two tissue slide images are unavailable for download from GDC Data Portal <!--DAT-1439-->
        * The raw and annotated VarScan VCF files for aliquot `TCGA-VR-A8ET-01A-11D-A403-09` are not available. These VCFs files will be replaced in a later release.<!--TT-602, DAT-1489-->
        * Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
        * Tumor grade property is not populated <!--SV-585-->
        * Progression_or_recurrence property is not populated <!--SV-584-->
    * TARGET projects
        * TARGET CGI BAMs in the Legacy Archive for the following aliquots should not be used because they were not repaired and concatenated into their original composite BAM files by CGHub.
            * TARGET-20-PASJGZ-04A-02D
            * TARGET-30-PAPTLY-01A-01D
            * TARGET-20-PAEIKD-09A-01D
            * TARGET-20-PASMYS-14A-02D
            * TARGET-20-PAMYAS-14A-02D
            * TARGET-10-PAPZST-09A-01D
        * 11 bam files for TARGET-NBL RNA-Seq are not available in the GDC Data portal <!--DAT-1476-->
        * There are 5051 TARGET files for which `experimental_strategy`, `data_format`, `platform`, and `data_subtype` are blank <!--SV-944-->
        * There are two cases with identical submitter_id `TARGET-10-PARUYU` <!--SV-940-->
        * Some TARGET cases are missing `days_to_last_follow_up` <!--SV-934-->
        * Some TARGET cases are missing `age_at_diagnosis` <!--SV-933-->
        * Some TARGET files are not connected to all related aliquots <!--SV-929-->
        * Samples of TARGET sample_type `Recurrent Blood Derived Cancer - Bone Marrow` are mislabeled as `Recurrent Blood Derived Cancer - Peripheral Blood`.  A workaround is to look at the sample barcode, which is -04 for `Recurrent Blood Derived Cancer - Bone Marrow`. (e.g. `TARGET-20-PAMYAS-04A-03R`) <!--SV-918-->
        * The latest TARGET data is not yet available at the GDC.  For the complete and latest data, please see the [TARGET Data Matrix](https://ocg.cancer.gov/programs/target/data-matrix).  Data that is not present or is not the most up to date includes:
            *  All microarray data and metadata
            *  All sequencing analyzed data and metadata
            *  1180 of 12063 sequencing runs of raw data
        * Demographic information for some TARGET patients is incorrect.  The correct information can be found in the associated clinical supplement file.  Impacted patients are TARGET-50-PAJNUS. <!--SV-710-->
        * No data from TARGET-MDLS is available.
    * Issues in the Legacy Archive
        * The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
        * SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
        * Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
        * SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
        * TARGET-MDLS cases do not have disease_type or primary_site populated <!--SV-939-->



## Data Release 27.0 Bug Fix

* __GDC Product__: Data
* __Release Date__: November 9, 2020

### New updates

1.  None, see bug fix section below.

A complete list of files for this release are listed for the GDC Data Portal and the GDC Legacy Archive are found below:

* [gdc_manifest_20201109_data_release_27.0_active.tsv.gz](gdc_manifest_20201109_data_release_27.0_active.tsv.gz)
* [gdc_manifest_20201109_data_release_27.0_legacy.tsv.gz](gdc_manifest_20201109_data_release_27.0_legacy.tsv.gz)

### Bugs Fixed Since Last Release

*  Some files in projects CGCI-BLGSP, CGCI-HTMCP-CC, and HCMI-CMDC were marked on the portal as controlled-access, when they were supposed to be open-access. These are now downloadable as open-access files.

### Known Issues and Workarounds

* The Copy Number Estimate files in GENIE are labeled on the portal as TXT while the files are actually in TSV format.  <!--DAT-2728-->
* Some tumor-only annotated VCFs (not raw VCFs) could have a small proportion of variants that appear twice.  Tumor-only annotated VCFs can be identified by searching for workflow "GATK4 MuTect2 Annotation" <!--SV-1425-->
* The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.
* Some miRNA files with QC failed reads were not swapped in DR11.0.  361 aliquots remain to be swapped in a later release <!--DAT-1589-->
* Mutation frequency may be underestimated when using MAF files for genes that overlap other genes.  This is because MAF files only record one gene per variant.
* Most intronic mutations are removed for MAF generation.  However, validated variants may rescue these in some cases.  Therefore intronic mutations in MAF files are not representative of those called by mutation callers.
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
* TCGA Projects
    * Incorrect information about treatment may be included for patients within TCGA-HNSC and TCGA-LGG.  Please refer to the clinical XML for accurate information on treatment <!--DAT-2264, DAT-2265-->
    * 74 Diagnostic TCGA slides are attached to a portion rather than a sample like the rest of the diagnostic slides. The reflects how these original samples were handled. <!--SV-1111-->
    * Two tissue slide images are unavailable for download from GDC Data Portal <!--DAT-1439-->
    * The raw and annotated VarScan VCF files for aliquot `TCGA-VR-A8ET-01A-11D-A403-09` are not available. These VCFs files will be replaced in a later release.<!--TT-602, DAT-1489-->
    * Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
    * Tumor grade property is not populated <!--SV-585-->
    * Progression_or_recurrence property is not populated <!--SV-584-->
* TARGET projects
    * TARGET CGI BAMs in the Legacy Archive for the following aliquots should not be used because they were not repaired and concatenated into their original composite BAM files by CGHub.
        * TARGET-20-PASJGZ-04A-02D
        * TARGET-30-PAPTLY-01A-01D
        * TARGET-20-PAEIKD-09A-01D
        * TARGET-20-PASMYS-14A-02D
        * TARGET-20-PAMYAS-14A-02D
        * TARGET-10-PAPZST-09A-01D
    * 11 bam files for TARGET-NBL RNA-Seq are not available in the GDC Data portal <!--DAT-1476-->
    * There are 5051 TARGET files for which `experimental_strategy`, `data_format`, `platform`, and `data_subtype` are blank <!--SV-944-->
    * There are two cases with identical submitter_id `TARGET-10-PARUYU` <!--SV-940-->
    * Some TARGET cases are missing `days_to_last_follow_up` <!--SV-934-->
    * Some TARGET cases are missing `age_at_diagnosis` <!--SV-933-->
    * Some TARGET files are not connected to all related aliquots <!--SV-929-->
    * Samples of TARGET sample_type `Recurrent Blood Derived Cancer - Bone Marrow` are mislabeled as `Recurrent Blood Derived Cancer - Peripheral Blood`.  A workaround is to look at the sample barcode, which is -04 for `Recurrent Blood Derived Cancer - Bone Marrow`. (e.g. `TARGET-20-PAMYAS-04A-03R`) <!--SV-918-->
    * The latest TARGET data is not yet available at the GDC.  For the complete and latest data, please see the [TARGET Data Matrix](https://ocg.cancer.gov/programs/target/data-matrix).  Data that is not present or is not the most up to date includes:
        *  All microarray data and metadata
        *  All sequencing analyzed data and metadata
        *  1180 of 12063 sequencing runs of raw data
    * Demographic information for some TARGET patients is incorrect.  The correct information can be found in the associated clinical supplement file.  Impacted patients are TARGET-50-PAJNUS. <!--SV-710-->
    * No data from TARGET-MDLS is available.
* Issues in the Legacy Archive
    * The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
    * SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
    * Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
    * SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
    * TARGET-MDLS cases do not have disease_type or primary_site populated <!--SV-939-->

## Data Release 27.0 <!--REQ-413-->

* __GDC Product__: Data
* __Release Date__: October 29, 2020

### New updates

1.  Initial release for the WGS variant calling pipeline. See the [documentation](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/#whole-genome-sequencing-variant-calling) on WGS variant calling for more details on the available files.  This includes data from the following projects:
    * CGCI-BLGSP <!--SPT7-75-->
    * CGCI-HTMCP-CC <!--SPT7-27-->
    * HCMI-CMDC <!--DAT-2948-->
2. RNA-Seq transcript fusion files are available for the following projects: <!--DAT-2947-->
    * CGCI-BLGSP
    * CGCI-HTMCP-CC
    * HCMI-CMDC
3. Aliquot level MAFs were released for CGCI-HTMCP-CC Targeted Sequencing variants. Open access MAFs are included. <!--SPT7-27-->
4. 17 new cases were released for the HCMI-CMDC project.  This includes WGS, WXS, and RNA-Seq data. <!--DAT-2949-->
5.  WGS alignments were released for 99 TCGA-LUAD cases (196 files). <!--DAT-1831-->
6. Therapeutic agents (treatment) and tumor stage (diagnosis) properties were migrated to remove deprecated values and better adhere to a standardized set of values. <!--DAT-2831--> <!--DAT-2880-->

A complete list of files for DR27.0 are listed for the GDC Data Portal and the GDC Legacy Archive are found below:

* [gdc_manifest_20201029_data_release_27.0_active.tsv.gz](gdc_manifest_20201029_data_release_27.0_active.tsv.gz)
* [gdc_manifest_20201029_data_release_27.0_legacy.tsv.gz](gdc_manifest_20201029_data_release_27.0_legacy.tsv.gz)

### Bugs Fixed Since Last Release

*  None

### Known Issues and Workarounds

* Some files in projects CGCI-BLGSP, CGCI-HTMCP-CC, and HCMI-CMDC are marked on the portal as controlled-access. These files are publicly downloadable using the Data Transfer Tool or API. All files from the following data types should be open-access within the previously specified projects: Biospecimen Supplement, Clinical Supplement, Gene Expression Quantification, Masked Somatic Mutation
* The Copy Number Estimate files in GENIE are labeled on the portal as TXT while the files are actually in TSV format.  <!--DAT-2728-->
* Some tumor-only annotated VCFs (not raw VCFs) could have a small proportion of variants that appear twice.  Tumor-only annotated VCFs can be identified by searching for workflow "GATK4 MuTect2 Annotation" <!--SV-1425-->
* The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.
* Some miRNA files with QC failed reads were not swapped in DR11.0.  361 aliquots remain to be swapped in a later release <!--DAT-1589-->
* Mutation frequency may be underestimated when using MAF files for genes that overlap other genes.  This is because MAF files only record one gene per variant.
* Most intronic mutations are removed for MAF generation.  However, validated variants may rescue these in some cases.  Therefore intronic mutations in MAF files are not representative of those called by mutation callers.
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
* TCGA Projects
    * Incorrect information about treatment may be included for patients within TCGA-HNSC and TCGA-LGG.  Please refer to the clinical XML for accurate information on treatment <!--DAT-2264, DAT-2265-->
    * 74 Diagnostic TCGA slides are attached to a portion rather than a sample like the rest of the diagnostic slides. The reflects how these original samples were handled. <!--SV-1111-->
    * Two tissue slide images are unavailable for download from GDC Data Portal <!--DAT-1439-->
    * The raw and annotated VarScan VCF files for aliquot `TCGA-VR-A8ET-01A-11D-A403-09` are not available. These VCFs files will be replaced in a later release.<!--TT-602, DAT-1489-->
    * Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
    * Tumor grade property is not populated <!--SV-585-->
    * Progression_or_recurrence property is not populated <!--SV-584-->
* TARGET projects
    * TARGET CGI BAMs in the Legacy Archive for the following aliquots should not be used because they were not repaired and concatenated into their original composite BAM files by CGHub.
        * TARGET-20-PASJGZ-04A-02D
        * TARGET-30-PAPTLY-01A-01D
        * TARGET-20-PAEIKD-09A-01D
        * TARGET-20-PASMYS-14A-02D
        * TARGET-20-PAMYAS-14A-02D
        * TARGET-10-PAPZST-09A-01D
    * 11 bam files for TARGET-NBL RNA-Seq are not available in the GDC Data portal <!--DAT-1476-->
    * There are 5051 TARGET files for which `experimental_strategy`, `data_format`, `platform`, and `data_subtype` are blank <!--SV-944-->
    * There are two cases with identical submitter_id `TARGET-10-PARUYU` <!--SV-940-->
    * Some TARGET cases are missing `days_to_last_follow_up` <!--SV-934-->
    * Some TARGET cases are missing `age_at_diagnosis` <!--SV-933-->
    * Some TARGET files are not connected to all related aliquots <!--SV-929-->
    * Samples of TARGET sample_type `Recurrent Blood Derived Cancer - Bone Marrow` are mislabeled as `Recurrent Blood Derived Cancer - Peripheral Blood`.  A workaround is to look at the sample barcode, which is -04 for `Recurrent Blood Derived Cancer - Bone Marrow`. (e.g. `TARGET-20-PAMYAS-04A-03R`) <!--SV-918-->
    * The latest TARGET data is not yet available at the GDC.  For the complete and latest data, please see the [TARGET Data Matrix](https://ocg.cancer.gov/programs/target/data-matrix).  Data that is not present or is not the most up to date includes:
        *  All microarray data and metadata
        *  All sequencing analyzed data and metadata
        *  1180 of 12063 sequencing runs of raw data
    * Demographic information for some TARGET patients is incorrect.  The correct information can be found in the associated clinical supplement file.  Impacted patients are TARGET-50-PAJNUS. <!--SV-710-->
    * No data from TARGET-MDLS is available.
* Issues in the Legacy Archive
    * The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
    * SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
    * Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
    * SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
    * TARGET-MDLS cases do not have disease_type or primary_site populated <!--SV-939-->

## Data Release 26.0 <!--REQ-409-->

* __GDC Product__: Data
* __Release Date__: September 8, 2020

### New updates

1.  New program released:
    * Count Me In (CMI)
        * CMI-ASC - The Angiosarcoma Project
            * RNA-Seq
            * WXS
        * CMI-MBC - The Metastatic Breast Cancer Project
            * RNA-Seq
            * WXS
2. Somatic mutations are now available on the exploration portal for the following projects:
    * MMRF-COMMPASS
    * TARGET-ALL-P3
    * TARGET-AML
    * TARGET-NBL
    * TARGET-WT
3. Primary sites and disease types were updated for multiple projects to correspond to GDC Dictionary updates.

A complete list of files for DR26.0 are listed for the GDC Data Portal and the GDC Legacy Archive are found below:

* [gdc_manifest_20200908_data_release_26.0_active.tsv.gz](gdc_manifest_20200908_data_release_26.0_active.tsv.gz)
* [gdc_manifest_20200908_data_release_26.0_legacy.tsv.gz](gdc_manifest_20200908_data_release_26.0_legacy.tsv.gz)

### Bugs Fixed Since Last Release

*  The CPTAC-3 head and neck cohort can now be queried by choosing the head and neck anatomic site on the GDC home page.

### Known Issues and Workarounds

* The Copy Number Estimate files in GENIE are labeled on the portal as TXT while the files are actually in TSV format.  <!--DAT-2728-->
* Some tumor-only annotated VCFs (not raw VCFs) could have a small proportion of variants that appear twice.  Tumor-only annotated VCFs can be identified by searching for workflow "GATK4 MuTect2 Annotation" <!--SV-1425-->
* The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.
* Some miRNA files with QC failed reads were not swapped in DR11.0.  361 aliquots remain to be swapped in a later release <!--DAT-1589-->
* Mutation frequency may be underestimated when using MAF files for genes that overlap other genes.  This is because MAF files only record one gene per variant.
* Most intronic mutations are removed for MAF generation.  However, validated variants may rescue these in some cases.  Therefore intronic mutations in MAF files are not representative of those called by mutation callers.
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
* TCGA Projects
    * Incorrect information about treatment may be included for patients within TCGA-HNSC and TCGA-LGG.  Please refer to the clinical XML for accurate information on treatment <!--DAT-2264, DAT-2265-->
    * 74 Diagnostic TCGA slides are attached to a portion rather than a sample like the rest of the diagnostic slides. The reflects how these original samples were handled. <!--SV-1111-->
    * Two tissue slide images are unavailable for download from GDC Data Portal <!--DAT-1439-->
    * The raw and annotated VarScan VCF files for aliquot `TCGA-VR-A8ET-01A-11D-A403-09` are not available. These VCFs files will be replaced in a later release.<!--TT-602, DAT-1489-->
    * Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
    * Tumor grade property is not populated <!--SV-585-->
    * Progression_or_recurrence property is not populated <!--SV-584-->
* TARGET projects
    * TARGET CGI BAMs in the Legacy Archive for the following aliquots should not be used because they were not repaired and concatenated into their original composite BAM files by CGHub.
        * TARGET-20-PASJGZ-04A-02D
        * TARGET-30-PAPTLY-01A-01D
        * TARGET-20-PAEIKD-09A-01D
        * TARGET-20-PASMYS-14A-02D
        * TARGET-20-PAMYAS-14A-02D
        * TARGET-10-PAPZST-09A-01D
    * 11 bam files for TARGET-NBL RNA-Seq are not available in the GDC Data portal <!--DAT-1476-->
    * There are 5051 TARGET files for which `experimental_strategy`, `data_format`, `platform`, and `data_subtype` are blank <!--SV-944-->
    * There are two cases with identical submitter_id `TARGET-10-PARUYU` <!--SV-940-->
    * Some TARGET cases are missing `days_to_last_follow_up` <!--SV-934-->
    * Some TARGET cases are missing `age_at_diagnosis` <!--SV-933-->
    * Some TARGET files are not connected to all related aliquots <!--SV-929-->
    * Samples of TARGET sample_type `Recurrent Blood Derived Cancer - Bone Marrow` are mislabeled as `Recurrent Blood Derived Cancer - Peripheral Blood`.  A workaround is to look at the sample barcode, which is -04 for `Recurrent Blood Derived Cancer - Bone Marrow`. (e.g. `TARGET-20-PAMYAS-04A-03R`) <!--SV-918-->
    * The latest TARGET data is not yet available at the GDC.  For the complete and latest data, please see the [TARGET Data Matrix](https://ocg.cancer.gov/programs/target/data-matrix).  Data that is not present or is not the most up to date includes:
        *  All microarray data and metadata
        *  All sequencing analyzed data and metadata
        *  1180 of 12063 sequencing runs of raw data
    * Demographic information for some TARGET patients is incorrect.  The correct information can be found in the associated clinical supplement file.  Impacted patients are TARGET-50-PAJNUS. <!--SV-710-->
    * No data from TARGET-MDLS is available.
* Issues in the Legacy Archive
    * The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
    * SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
    * Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
    * SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
    * TARGET-MDLS cases do not have disease_type or primary_site populated <!--SV-939-->

    ## Data Release 25.0 <!--REQ-403-->

    * __GDC Product__: Data
    * __Release Date__: July 22, 2020

    ### New updates

    1.  New data types released:
        * RNA-Seq Transcript Fusion files were released for the following projects:
            * TARGET-ALL-P1
            * TARGET-ALL-P2
            * TARGET-ALL-P3
            * TARGET-CCSK
            * TARGET-NBL
            * TARGET-OS
            * TARGET-RT
            * TARGET-WT
        * The msi_status and msi_score properties can be queried on the GDC Portal for the CPTAC-3 project.  
            * To query for these fields: go to the [GDC Repository](https://portal.gdc.cancer.gov/repository), click on "Add a File Filter" at the top left of the screen, type msi_score or msi_status in the field, and click on "msi_score" or "msi_status".  This should bring up the corresponding filters to use on the portal.
    2. 108 cases from the CPTAC-3 LSCC Cohort were released. Includes the following data types:
        * WXS
        * WGS
        * RNA-Seq
        * miRNA-Seq
    3. Aliquot level MAFs were released for MMRF-COMMPASS WXS variants. Open access MAFs are included.
    4. HCMI-CMDC open-access somatic mutations were released to the [Exploration Portal](https://portal.gdc.cancer.gov/exploration).  

    A complete list of files for DR25.0 are listed for the GDC Data Portal and the GDC Legacy Archive are found below:

    * [gdc_manifest_20200722_data_release_25.0_active.tsv.gz](gdc_manifest_20200722_data_release_25.0_active.tsv.gz)
    * [gdc_manifest_20200722_data_release_25.0_legacy.tsv.gz](gdc_manifest_20200722_data_release_25.0_legacy.tsv.gz)

    ### Bugs Fixed Since Last Release

    *  A few supplements from CGCI-BLGSP are now associated with their correct versions. <!--SV-1709-->

    ### Known Issues and Workarounds

    * Currently the CPTAC-3 HNSCC cohort does not appear when the "Head and Neck" primary site is selected from the GDC home page. This cohort can be queried by clicking [here](https://portal.gdc.cancer.gov/repository?facetTab=cases&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.primary_site%22%2C%22value%22%3A%5B%22other%20and%20ill-defined%20sites%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%5B%22CPTAC-3%22%5D%7D%7D%5D%7D)
    * The Copy Number Estimate files in GENIE are labeled on the portal as TXT while the files are actually in TSV format.  <!--DAT-2728-->
    * Some tumor-only annotated VCFs (not raw VCFs) could have a small proportion of variants that appear twice.  Tumor-only annotated VCFs can be identified by searching for workflow "GATK4 MuTect2 Annotation" <!--SV-1425-->
    * The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.
    * Some miRNA files with QC failed reads were not swapped in DR11.0.  361 aliquots remain to be swapped in a later release <!--DAT-1589-->
    * Mutation frequency may be underestimated when using MAF files for genes that overlap other genes.  This is because MAF files only record one gene per variant.
    * Most intronic mutations are removed for MAF generation.  However, validated variants may rescue these in some cases.  Therefore intronic mutations in MAF files are not representative of those called by mutation callers.
    * Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
    * BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
    * Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
    * TCGA Projects
        * Incorrect information about treatment may be included for patients within TCGA-HNSC and TCGA-LGG.  Please refer to the clinical XML for accurate information on treatment <!--DAT-2264, DAT-2265-->
        * 74 Diagnostic TCGA slides are attached to a portion rather than a sample like the rest of the diagnostic slides. The reflects how these original samples were handled. <!--SV-1111-->
        * Two tissue slide images are unavailable for download from GDC Data Portal <!--DAT-1439-->
        * The raw and annotated VarScan VCF files for aliquot `TCGA-VR-A8ET-01A-11D-A403-09` are not available. These VCFs files will be replaced in a later release.<!--TT-602, DAT-1489-->
        * Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
        * Tumor grade property is not populated <!--SV-585-->
        * Progression_or_recurrence property is not populated <!--SV-584-->
    * TARGET projects
        * TARGET CGI BAMs in the Legacy Archive for the following aliquots should not be used because they were not repaired and concatenated into their original composite BAM files by CGHub.
            * TARGET-20-PASJGZ-04A-02D
            * TARGET-30-PAPTLY-01A-01D
            * TARGET-20-PAEIKD-09A-01D
            * TARGET-20-PASMYS-14A-02D
            * TARGET-20-PAMYAS-14A-02D
            * TARGET-10-PAPZST-09A-01D
        * 11 bam files for TARGET-NBL RNA-Seq are not available in the GDC Data portal <!--DAT-1476-->
        * There are 5051 TARGET files for which `experimental_strategy`, `data_format`, `platform`, and `data_subtype` are blank <!--SV-944-->
        * There are two cases with identical submitter_id `TARGET-10-PARUYU` <!--SV-940-->
        * Some TARGET cases are missing `days_to_last_follow_up` <!--SV-934-->
        * Some TARGET cases are missing `age_at_diagnosis` <!--SV-933-->
        * Some TARGET files are not connected to all related aliquots <!--SV-929-->
        * Samples of TARGET sample_type `Recurrent Blood Derived Cancer - Bone Marrow` are mislabeled as `Recurrent Blood Derived Cancer - Peripheral Blood`.  A workaround is to look at the sample barcode, which is -04 for `Recurrent Blood Derived Cancer - Bone Marrow`. (e.g. `TARGET-20-PAMYAS-04A-03R`) <!--SV-918-->
        * The latest TARGET data is not yet available at the GDC.  For the complete and latest data, please see the [TARGET Data Matrix](https://ocg.cancer.gov/programs/target/data-matrix).  Data that is not present or is not the most up to date includes:
            *  All microarray data and metadata
            *  All sequencing analyzed data and metadata
            *  1180 of 12063 sequencing runs of raw data
        * Demographic information for some TARGET patients is incorrect.  The correct information can be found in the associated clinical supplement file.  Impacted patients are TARGET-50-PAJNUS. <!--SV-710-->
        * No data from TARGET-MDLS is available.
    * Issues in the Legacy Archive
        * The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
        * SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
        * Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
        * SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
        * TARGET-MDLS cases do not have disease_type or primary_site populated <!--SV-939-->


## Data Release 24.0 <!--REQ-401-->

* __GDC Product__: Data
* __Release Date__: May 7, 2020

### New updates

1.  New project released: CGCI-HTMCP-CC - HIV+ Tumor Molecular Characterization Project - Cervical Cancer
    * RNA-Seq: Alignments and gene expression levels
    * miRNA-Seq: Alignments and miRNA expression levels
    * WGS: Alignments
    * Targeted Sequencing: Alignments

2. 110 new cases were released from the HNSCC cohort of CPTAC-3.  This includes WXS, WGS, RNA-Seq and miRNA-Seq data.

3. Aliquot-level WXS MAFs are now available from the following projects:
    * CPTAC-2
    * CPTAC-3

A complete list of files for DR24.0 are listed for the GDC Data Portal and the GDC Legacy Archive are found below:

* [gdc_manifest_20200507_data_release_24.0_active.tsv.gz](gdc_manifest_20200507_data_release_24.0_active.tsv.gz)
* [gdc_manifest_20200507_data_release_24.0_legacy.tsv.gz](gdc_manifest_20200507_data_release_24.0_legacy.tsv.gz)


### Bugs Fixed Since Last Release

*  None

### Known Issues and Workarounds

* Currently the CPTAC-3 HNSCC cohort does not appear when the "Head and Neck" primary site is selected from the GDC home page. This cohort can be queried by clicking [here](https://portal.gdc.cancer.gov/repository?facetTab=cases&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.primary_site%22%2C%22value%22%3A%5B%22other%20and%20ill-defined%20sites%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%5B%22CPTAC-3%22%5D%7D%7D%5D%7D)
* The Copy Number Estimate files in GENIE are labeled on the portal as TXT while the files are actually in TSV format.  <!--DAT-2728-->
* Some tumor-only annotated VCFs (not raw VCFs) could have a small proportion of variants that appear twice.  Tumor-only annotated VCFs can be identified by searching for workflow "GATK4 MuTect2 Annotation" <!--SV-1425-->
* The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.
* Some miRNA files with QC failed reads were not swapped in DR11.0.  361 aliquots remain to be swapped in a later release <!--DAT-1589-->
* Mutation frequency may be underestimated when using MAF files for genes that overlap other genes.  This is because MAF files only record one gene per variant.
* Most intronic mutations are removed for MAF generation.  However, validated variants may rescue these in some cases.  Therefore intronic mutations in MAF files are not representative of those called by mutation callers.
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
* TCGA Projects
    * Incorrect information about treatment may be included for patients within TCGA-HNSC and TCGA-LGG.  Please refer to the clinical XML for accurate information on treatment <!--DAT-2264, DAT-2265-->
    * 74 Diagnostic TCGA slides are attached to a portion rather than a sample like the rest of the diagnostic slides. The reflects how these original samples were handled. <!--SV-1111-->
    * Two tissue slide images are unavailable for download from GDC Data Portal <!--DAT-1439-->
    * The raw and annotated VarScan VCF files for aliquot `TCGA-VR-A8ET-01A-11D-A403-09` are not available. These VCFs files will be replaced in a later release.<!--TT-602, DAT-1489-->
    * Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
    * Tumor grade property is not populated <!--SV-585-->
    * Progression_or_recurrence property is not populated <!--SV-584-->
* TARGET projects
    * TARGET CGI BAMs in the Legacy Archive for the following aliquots should not be used because they were not repaired and concatenated into their original composite BAM files by CGHub.
        * TARGET-20-PASJGZ-04A-02D
        * TARGET-30-PAPTLY-01A-01D
        * TARGET-20-PAEIKD-09A-01D
        * TARGET-20-PASMYS-14A-02D
        * TARGET-20-PAMYAS-14A-02D
        * TARGET-10-PAPZST-09A-01D
    * 11 bam files for TARGET-NBL RNA-Seq are not available in the GDC Data portal <!--DAT-1476-->
    * There are 5051 TARGET files for which `experimental_strategy`, `data_format`, `platform`, and `data_subtype` are blank <!--SV-944-->
    * There are two cases with identical submitter_id `TARGET-10-PARUYU` <!--SV-940-->
    * Some TARGET cases are missing `days_to_last_follow_up` <!--SV-934-->
    * Some TARGET cases are missing `age_at_diagnosis` <!--SV-933-->
    * Some TARGET files are not connected to all related aliquots <!--SV-929-->
    * Samples of TARGET sample_type `Recurrent Blood Derived Cancer - Bone Marrow` are mislabeled as `Recurrent Blood Derived Cancer - Peripheral Blood`.  A workaround is to look at the sample barcode, which is -04 for `Recurrent Blood Derived Cancer - Bone Marrow`. (e.g. `TARGET-20-PAMYAS-04A-03R`) <!--SV-918-->
    * The latest TARGET data is not yet available at the GDC.  For the complete and latest data, please see the [TARGET Data Matrix](https://ocg.cancer.gov/programs/target/data-matrix).  Data that is not present or is not the most up to date includes:
        *  All microarray data and metadata
        *  All sequencing analyzed data and metadata
        *  1180 of 12063 sequencing runs of raw data
    * Demographic information for some TARGET patients is incorrect.  The correct information can be found in the associated clinical supplement file.  Impacted patients are TARGET-50-PAJNUS. <!--SV-710-->
    * No data from TARGET-MDLS is available.
* Issues in the Legacy Archive
    * The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
    * SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
    * Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
    * SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
    * TARGET-MDLS cases do not have disease_type or primary_site populated <!--SV-939-->

## Data Release 23.0 <!--REQ-397-->

* __GDC Product__: Data
* __Release Date__: April 7, 2020

### New updates

1.  New data types released:
    * Aliquot-level MAFs: MAF Files with mutations derived from one tumor/normal pair <!--DAT-2732-->
        * HCMI-CMDC
        * TARGET-ALL-P2
        * TARGET-ALL-P3
        * TARGET-AML
        * TARGET-NBL
        * TARGET-OS
        * TARGET-WT
        * Note: Previously released TARGET project level MAFs can be downloaded with the following manifest: [TARGET_Project-Level-MAF_GDC-Manifest.txt](TARGET_Project-Level-MAF_GDC-Manifest.txt)
    * Copy number segment and estimate files from SNP6 ASCAT <!--BINF-293-->
        * All TCGA Projects
        * TARGET-ALL-P2
        * TARGET-AML

2.  To accommodate users who prefer to use project-level MAFs, a MAF aggregation tool was developed by the GDC: <!--DEV-89-->
    * [Github Release](https://github.com/NCI-GDC/gdc-maf-tool/releases)

3.  New RNA-Seq data was released from HCMI-CMDC for nine additional cases.  

4.  Clinical updates were performed for the following projects
    * CGCI-BLGSP
    * HCMI-CMDC <!--SV-1660-->
    * WCDT-MCRPC <!--DAT-2782-->

A complete list of files for DR23.0 are listed for the GDC Data Portal and the GDC Legacy Archive are found below:

* [gdc_manifest_20200407_data_release_23.0_active.tsv.gz](gdc_manifest_20200407_data_release_23.0_active.tsv.gz)
* [gdc_manifest_20200407_data_release_23.0_legacy.tsv.gz](gdc_manifest_20200407_data_release_23.0_legacy.tsv.gz)


### Bugs Fixed Since Last Release

*  The 6 HCMI-CMDC cases without clinical data now have clinical data. <!--SV-1660-->
*  Most of the "associated_entities" fields in CGCI-BLGSP were not populated correct, this has been resolved.  <!--SV-1701-->

### Known Issues and Workarounds

* The Copy Number Estimate files in GENIE are labeled on the portal as TXT while the files are actually in TSV format.  <!--DAT-2728-->
* Some tumor-only annotated VCFs (not raw VCFs) could have a small proportion of variants that appear twice.  Tumor-only annotated VCFs can be identified by searching for workflow "GATK4 MuTect2 Annotation" <!--SV-1425-->
* The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.
* Some miRNA files with QC failed reads were not swapped in DR11.0.  361 aliquots remain to be swapped in a later release <!--DAT-1589-->
* Mutation frequency may be underestimated when using MAF files for genes that overlap other genes.  This is because MAF files only record one gene per variant.
* Most intronic mutations are removed for MAF generation.  However, validated variants may rescue these in some cases.  Therefore intronic mutations in MAF files are not representative of those called by mutation callers.
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
* TCGA Projects
    * Incorrect information about treatment may be included for patients within TCGA-HNSC and TCGA-LGG.  Please refer to the clinical XML for accurate information on treatment <!--DAT-2264, DAT-2265-->
    * 74 Diagnostic TCGA slides are attached to a portion rather than a sample like the rest of the diagnostic slides. The reflects how these original samples were handled. <!--SV-1111-->
    * Two tissue slide images are unavailable for download from GDC Data Portal <!--DAT-1439-->
    * The raw and annotated VarScan VCF files for aliquot `TCGA-VR-A8ET-01A-11D-A403-09` are not available. These VCFs files will be replaced in a later release.<!--TT-602, DAT-1489-->
    * Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
    * Tumor grade property is not populated <!--SV-585-->
    * Progression_or_recurrence property is not populated <!--SV-584-->
* TARGET projects
    * TARGET CGI BAMs in the Legacy Archive for the following aliquots should not be used because they were not repaired and concatenated into their original composite BAM files by CGHub.
        * TARGET-20-PASJGZ-04A-02D
        * TARGET-30-PAPTLY-01A-01D
        * TARGET-20-PAEIKD-09A-01D
        * TARGET-20-PASMYS-14A-02D
        * TARGET-20-PAMYAS-14A-02D
        * TARGET-10-PAPZST-09A-01D
    * 11 bam files for TARGET-NBL RNA-Seq are not available in the GDC Data portal <!--DAT-1476-->
    * There are 5051 TARGET files for which `experimental_strategy`, `data_format`, `platform`, and `data_subtype` are blank <!--SV-944-->
    * There are two cases with identical submitter_id `TARGET-10-PARUYU` <!--SV-940-->
    * Some TARGET cases are missing `days_to_last_follow_up` <!--SV-934-->
    * Some TARGET cases are missing `age_at_diagnosis` <!--SV-933-->
    * Some TARGET files are not connected to all related aliquots <!--SV-929-->
    * Samples of TARGET sample_type `Recurrent Blood Derived Cancer - Bone Marrow` are mislabeled as `Recurrent Blood Derived Cancer - Peripheral Blood`.  A workaround is to look at the sample barcode, which is -04 for `Recurrent Blood Derived Cancer - Bone Marrow`. (e.g. `TARGET-20-PAMYAS-04A-03R`) <!--SV-918-->
    * The latest TARGET data is not yet available at the GDC.  For the complete and latest data, please see the [TARGET Data Matrix](https://ocg.cancer.gov/programs/target/data-matrix).  Data that is not present or is not the most up to date includes:
        *  All microarray data and metadata
        *  All sequencing analyzed data and metadata
        *  1180 of 12063 sequencing runs of raw data
    * Demographic information for some TARGET patients is incorrect.  The correct information can be found in the associated clinical supplement file.  Impacted patients are TARGET-50-PAJNUS. <!--SV-710-->
    * No data from TARGET-MDLS is available.
* Issues in the Legacy Archive
    * The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
    * SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
    * Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
    * SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
    * TARGET-MDLS cases do not have disease_type or primary_site populated <!--SV-939-->

## Data Release 22.0 <!--REQ-396-->

* __GDC Product__: Data
* __Release Date__: January 16, 2020

### New updates

1.  New projects released:
    * WCDT-MCRPC - Genomic Characterization of Metastatic Castration Resistant Prostate Cancer (phs001648)
        * RNA-Seq; WGS Data <!--SPT7-4--> <!--SPT7-38-->
2. New data from HCMI-CMDC <!--DAT-2758--> <!--SPT7-33-->
    * 16 New Cases
    * Includes WXS, WGS, and RNA-Seq data
3. New data from CPTAC-3 <!--DAT-2723-->
    * 108 New Cases
    * Includes WXS, WGS, and RNA-Seq data
    * miRNA-Seq data for currently released cases

A complete list of files for DR22.0 are listed for the GDC Data Portal and the GDC Legacy Archive are found below:

* [gdc_manifest_20200116_data_release_22.0_active.tsv.gz](gdc_manifest_20200116_data_release_22.0_active.tsv.gz)
* [gdc_manifest_20200116_data_release_22.0_legacy.tsv.gz](gdc_manifest_20200116_data_release_22.0_legacy.tsv.gz)


### Bugs Fixed Since Last Release

*  None

### Known Issues and Workarounds

* The Copy Number Estimate files in GENIE are labeled on the portal as TXT while the files are actually in TSV format.  <!--DAT-2728-->
* 6 of the HCMI-CMDC cases are missing clinical nodes <!--SV-1660-->
    * HCM-CSHL-0060-C18
    * HCM-CSHL-0089-C25
    * HCM-CSHL-0090-C25
    * HCM-CSHL-0092-C25
    * HCM-CSHL-0091-C25
    * HCM-CSHL-0057-C18
* Some tumor-only annotated VCFs (not raw VCFs) could have a small proportion of variants that appear twice.  Tumor-only annotated VCFs can be identified by searching for workflow "GATK4 MuTect2 Annotation" <!--SV-1425-->
* The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.
* Some miRNA files with QC failed reads were not swapped in DR11.0.  361 aliquots remain to be swapped in a later release <!--DAT-1589-->
* Mutation frequency may be underestimated when using MAF files for genes that overlap other genes.  This is because MAF files only record one gene per variant.
* Most intronic mutations are removed for MAF generation.  However, validated variants may rescue these in some cases.  Therefore intronic mutations in MAF files are not representative of those called by mutation callers.
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
* TCGA Projects
    * Incorrect information about treatment may be included for patients within TCGA-HNSC and TCGA-LGG.  Please refer to the clinical XML for accurate information on treatment <!--DAT-2264, DAT-2265-->
    * 74 Diagnostic TCGA slides are attached to a portion rather than a sample like the rest of the diagnostic slides. The reflects how these original samples were handled. <!--SV-1111-->
    * Two tissue slide images are unavailable for download from GDC Data Portal <!--DAT-1439-->
    * The raw and annotated VarScan VCF files for aliquot `TCGA-VR-A8ET-01A-11D-A403-09` are not available. These VCFs files will be replaced in a later release.<!--TT-602, DAT-1489-->
    * Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
    * Tumor grade property is not populated <!--SV-585-->
    * Progression_or_recurrence property is not populated <!--SV-584-->
* TARGET projects
    * TARGET CGI BAMs in the Legacy Archive for the following aliquots should not be used because they were not repaired and concatenated into their original composite BAM files by CGHub.
        * TARGET-20-PASJGZ-04A-02D
        * TARGET-30-PAPTLY-01A-01D
        * TARGET-20-PAEIKD-09A-01D
        * TARGET-20-PASMYS-14A-02D
        * TARGET-20-PAMYAS-14A-02D
        * TARGET-10-PAPZST-09A-01D
    * 11 bam files for TARGET-NBL RNA-Seq are not available in the GDC Data portal <!--DAT-1476-->
    * There are 5051 TARGET files for which `experimental_strategy`, `data_format`, `platform`, and `data_subtype` are blank <!--SV-944-->
    * There are two cases with identical submitter_id `TARGET-10-PARUYU` <!--SV-940-->
    * Some TARGET cases are missing `days_to_last_follow_up` <!--SV-934-->
    * Some TARGET cases are missing `age_at_diagnosis` <!--SV-933-->
    * Some TARGET files are not connected to all related aliquots <!--SV-929-->
    * Samples of TARGET sample_type `Recurrent Blood Derived Cancer - Bone Marrow` are mislabeled as `Recurrent Blood Derived Cancer - Peripheral Blood`.  A workaround is to look at the sample barcode, which is -04 for `Recurrent Blood Derived Cancer - Bone Marrow`. (e.g. `TARGET-20-PAMYAS-04A-03R`) <!--SV-918-->
    * The latest TARGET data is not yet available at the GDC.  For the complete and latest data, please see the [TARGET Data Matrix](https://ocg.cancer.gov/programs/target/data-matrix).  Data that is not present or is not the most up to date includes:
        *  All microarray data and metadata
        *  All sequencing analyzed data and metadata
        *  1180 of 12063 sequencing runs of raw data
    * Demographic information for some TARGET patients is incorrect.  The correct information can be found in the associated clinical supplement file.  Impacted patients are TARGET-50-PAJNUS. <!--SV-710-->
    * No data from TARGET-MDLS is available.
* Issues in the Legacy Archive
    * The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
    * SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
    * Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
    * SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
    * TARGET-MDLS cases do not have disease_type or primary_site populated <!--SV-939-->

## Data Release 21.0 <!--REQ-396-->

* __GDC Product__: Data
* __Release Date__: December 10, 2019

### New updates

1.  New projects released:
    *  GENIE - AACR Project Genomics Evidence Neoplasia Information Exchange (phs001337) <!--SPT-173-->
        *  Includes Targeted Sequencing, Transcript Fusion, Copy Number Estimate from GENIE 5.0
    *  AACR Project GENIE is divided by sequencing center:
        * GENIE-MSK
	    * GENIE-DFCI
	    * GENIE-MDA
	    * GENIE-JHU
	    * GENIE-UHN
	    * GENIE-VICC
	    * GENIE-GRCC
	    * GENIE-NKI

A complete list of files for DR21.0 are listed for the GDC Data Portal and the GDC Legacy Archive are found below:

* [gdc_manifest_20191210_data_release_21.0_active.txt.gz](gdc_manifest_20191210_data_release_21.0_active.txt.gz)
* [gdc_manifest_20191210_data_release_21.0_legacy.txt.gz](gdc_manifest_20191210_data_release_21.0_legacy.txt.gz)


### Bugs Fixed Since Last Release

*  None

### Known Issues and Workarounds

* The Copy Number Estimate files in GENIE are labeled on the portal as TXT while the files are actually in TSV format.  <!--DAT-2728-->
* Some tumor-only annotated VCFs (not raw VCFs) could have a small proportion of variants that appear twice.  Tumor-only annotated VCFs can be identified by searching for workflow "GATK4 MuTect2 Annotation" <!--SV-1425-->
* The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.
* Some miRNA files with QC failed reads were not swapped in DR11.0.  361 aliquots remain to be swapped in a later release <!--DAT-1589-->
* Mutation frequency may be underestimated when using MAF files for genes that overlap other genes.  This is because MAF files only record one gene per variant.
* Most intronic mutations are removed for MAF generation.  However, validated variants may rescue these in some cases.  Therefore intronic mutations in MAF files are not representative of those called by mutation callers.
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
* TCGA Projects
    * Incorrect information about treatment may be included for patients within TCGA-HNSC and TCGA-LGG.  Please refer to the clinical XML for accurate information on treatment <!--DAT-2264, DAT-2265-->
    * 74 Diagnostic TCGA slides are attached to a portion rather than a sample like the rest of the diagnostic slides. The reflects how these original samples were handled. <!--SV-1111-->
    * Two tissue slide images are unavailable for download from GDC Data Portal <!--DAT-1439-->
    * The raw and annotated VarScan VCF files for aliquot `TCGA-VR-A8ET-01A-11D-A403-09` are not available. These VCFs files will be replaced in a later release.<!--TT-602, DAT-1489-->
    * Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
    * Tumor grade property is not populated <!--SV-585-->
    * Progression_or_recurrence property is not populated <!--SV-584-->
* TARGET projects
    * TARGET CGI BAMs in the Legacy Archive for the following aliquots should not be used because they were not repaired and concatenated into their original composite BAM files by CGHub.
        * TARGET-20-PASJGZ-04A-02D
        * TARGET-30-PAPTLY-01A-01D
        * TARGET-20-PAEIKD-09A-01D
        * TARGET-20-PASMYS-14A-02D
        * TARGET-20-PAMYAS-14A-02D
        * TARGET-10-PAPZST-09A-01D
    * 11 bam files for TARGET-NBL RNA-Seq are not available in the GDC Data portal <!--DAT-1476-->
    * There are 5051 TARGET files for which `experimental_strategy`, `data_format`, `platform`, and `data_subtype` are blank <!--SV-944-->
    * There are two cases with identical submitter_id `TARGET-10-PARUYU` <!--SV-940-->
    * Some TARGET cases are missing `days_to_last_follow_up` <!--SV-934-->
    * Some TARGET cases are missing `age_at_diagnosis` <!--SV-933-->
    * Some TARGET files are not connected to all related aliquots <!--SV-929-->
    * Samples of TARGET sample_type `Recurrent Blood Derived Cancer - Bone Marrow` are mislabeled as `Recurrent Blood Derived Cancer - Peripheral Blood`.  A workaround is to look at the sample barcode, which is -04 for `Recurrent Blood Derived Cancer - Bone Marrow`. (e.g. `TARGET-20-PAMYAS-04A-03R`) <!--SV-918-->
    * The latest TARGET data is not yet available at the GDC.  For the complete and latest data, please see the [TARGET Data Matrix](https://ocg.cancer.gov/programs/target/data-matrix).  Data that is not present or is not the most up to date includes:
        *  All microarray data and metadata
        *  All sequencing analyzed data and metadata
        *  1180 of 12063 sequencing runs of raw data
    * Demographic information for some TARGET patients is incorrect.  The correct information can be found in the associated clinical supplement file.  Impacted patients are TARGET-50-PAJNUS. <!--SV-710-->
    * No data from TARGET-MDLS is available.
* Issues in the Legacy Archive
    * The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
    * SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
    * Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
    * SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
    * TARGET-MDLS cases do not have disease_type or primary_site populated <!--SV-939-->

## Data Release 20.0 <!--REQ-392-->

* __GDC Product__: Data
* __Release Date__: November 11, 2019

### New updates

1.  New projects released:
    *  CPTAC-2 - CPTAC Proteogenomic Confirmatory Study (phs000892) <!--SPT-825-->
        *  Includes WXS, RNA-Seq, and miRNA-Seq
    *  OHSU-CNL - Genomic landscape of Neutrophilic Leukemias of Ambiguous Diagnosis (phs001799) <!--SPT-755-->
        *  Includes WXS and RNA-Seq
        *  No VCF files will be included at this time.  They will follow in a later release.
2.  New TARGET data released
    *  TARGET-OS: WGS, WXS <!--SPT-394, SPT-226-->
    *  TARGET-NBL: WGS <!--SPT-412-->
    *  TARGET-AML: miRNA <!--SV-1518-->
3.  CGCI-BLGSP miRNA-Seq released <!-- SPT-441-->


A complete list of files for DR20.0 are listed for the GDC Data Portal and the GDC Legacy Archive are found below:


* [gdc_manifest_20191111_data_release_20.0_active.txt.gz](gdc_manifest_20191111_data_release_20.0_active.txt.gz)
* [gdc_manifest_20191111_data_release_20.0_legacy.txt.gz](gdc_manifest_20191111_data_release_20.0_legacy.txt.gz)


### Bugs Fixed Since Last Release

*  None

### Known Issues and Workarounds

* Some tumor-only annotated VCFs (not raw VCFs) could have a small proportion of variants that appear twice.  Tumor-only annotated VCFs can be identified by searching for workflow "GATK4 MuTect2 Annotation" <!--SV-1425-->
* The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.
* Some miRNA files with QC failed reads were not swapped in DR11.0.  361 aliquots remain to be swapped in a later release <!--DAT-1589-->
* Mutation frequency may be underestimated when using MAF files for genes that overlap other genes.  This is because MAF files only record one gene per variant.
* Most intronic mutations are removed for MAF generation.  However, validated variants may rescue these in some cases.  Therefore intronic mutations in MAF files are not representative of those called by mutation callers.
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
* TCGA Projects
    * Incorrect information about treatment may be included for patients within TCGA-HNSC and TCGA-LGG.  Please refer to the clinical XML for accurate information on treatment <!--DAT-2264, DAT-2265-->
    * 74 Diagnostic TCGA slides are attached to a portion rather than a sample like the rest of the diagnostic slides. The reflects how these original samples were handled. <!--SV-1111-->
    * Two tissue slide images are unavailable for download from GDC Data Portal <!--DAT-1439-->
    * The raw and annotated VarScan VCF files for aliquot `TCGA-VR-A8ET-01A-11D-A403-09` are not available. These VCFs files will be replaced in a later release.<!--TT-602, DAT-1489-->
    * Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
    * Tumor grade property is not populated <!--SV-585-->
    * Progression_or_recurrence property is not populated <!--SV-584-->
* TARGET projects
    * TARGET CGI BAMs in the Legacy Archive for the following aliquots should not be used because they were not repaired and concatenated into their original composite BAM files by CGHub.
        * TARGET-20-PASJGZ-04A-02D
        * TARGET-30-PAPTLY-01A-01D
        * TARGET-20-PAEIKD-09A-01D
        * TARGET-20-PASMYS-14A-02D
        * TARGET-20-PAMYAS-14A-02D
        * TARGET-10-PAPZST-09A-01D
    * 11 bam files for TARGET-NBL RNA-Seq are not available in the GDC Data portal <!--DAT-1476-->
    * There are 5051 TARGET files for which `experimental_strategy`, `data_format`, `platform`, and `data_subtype` are blank <!--SV-944-->
    * There are two cases with identical submitter_id `TARGET-10-PARUYU` <!--SV-940-->
    * Some TARGET cases are missing `days_to_last_follow_up` <!--SV-934-->
    * Some TARGET cases are missing `age_at_diagnosis` <!--SV-933-->
    * Some TARGET files are not connected to all related aliquots <!--SV-929-->
    * Samples of TARGET sample_type `Recurrent Blood Derived Cancer - Bone Marrow` are mislabeled as `Recurrent Blood Derived Cancer - Peripheral Blood`.  A workaround is to look at the sample barcode, which is -04 for `Recurrent Blood Derived Cancer - Bone Marrow`. (e.g. `TARGET-20-PAMYAS-04A-03R`) <!--SV-918-->
    * The latest TARGET data is not yet available at the GDC.  For the complete and latest data, please see the [TARGET Data Matrix](https://ocg.cancer.gov/programs/target/data-matrix).  Data that is not present or is not the most up to date includes:
        *  All microarray data and metadata
        *  All sequencing analyzed data and metadata
        *  1180 of 12063 sequencing runs of raw data
    * Demographic information for some TARGET patients is incorrect.  The correct information can be found in the associated clinical supplement file.  Impacted patients are TARGET-50-PAJNUS. <!--SV-710-->
    * No data from TARGET-MDLS is available.
* Issues in the Legacy Archive
    * The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
    * SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
    * Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
    * SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
    * TARGET-MDLS cases do not have disease_type or primary_site populated <!--SV-939-->




## Data Release 19.1

* __GDC Product__: Data
* __Release Date__: November 6, 2019

### New updates

* The following cases are no longer available in the GDC Data Portal.  They had no data files associated with them in DR 19 so there are no changes in file availability in this release.
    * TARGET-00-NAAENF
    * TARGET-00-NAAENG
    * TARGET-00-NAAENH
    * TARGET-00-NAAENI
    * TARGET-00-NAAENJ
    * TARGET-00-NAAENK
    * TARGET-00-NAAENL
    * TARGET-00-NAAENM
    * TARGET-00-NAAENN
    * TARGET-00-NAAENP
    * TARGET-00-NAAENR
    * TARGET-00-NAAEPE

A complete list of files for DR19.1 are listed for the GDC Data Portal and the GDC Legacy Archive are found below:

* [gdc_manifest_20190917_data_release_19.0_active.txt.gz](gdc_manifest_20190917_data_release_19.0_active.txt.gz)
* [gdc_manifest_20190917_data_release_19.0_legacy.txt.gz](gdc_manifest_20190917_data_release_19.0_legacy.txt.gz)


### Bugs Fixed Since Last Release

*  None

### Known Issues and Workarounds

* Some tumor-only annotated VCFs (not raw VCFs) could have a small proportion of variants that appear twice.  Tumor-only annotated VCFs can be identified by searching for workflow "GATK4 MuTect2 Annotation" <!--SV-1425-->
* The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.
* Some miRNA files with QC failed reads were not swapped in DR11.0.  361 aliquots remain to be swapped in a later release <!--DAT-1589-->
* Mutation frequency may be underestimated when using MAF files for genes that overlap other genes.  This is because MAF files only record one gene per variant.
* Most intronic mutations are removed for MAF generation.  However, validated variants may rescue these in some cases.  Therefore intronic mutations in MAF files are not representative of those called by mutation callers.
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
* TCGA Projects
    * Incorrect information about treatment may be included for patients within TCGA-HNSC and TCGA-LGG.  Please refer to the clinical XML for accurate information on treatment <!--DAT-2264, DAT-2265-->
    * 74 Diagnostic TCGA slides are attached to a portion rather than a sample like the rest of the diagnostic slides. The reflects how these original samples were handled. <!--SV-1111-->
    * Two tissue slide images are unavailable for download from GDC Data Portal <!--DAT-1439-->
    * The raw and annotated VarScan VCF files for aliquot `TCGA-VR-A8ET-01A-11D-A403-09` are not available. These VCFs files will be replaced in a later release.<!--TT-602, DAT-1489-->
    * Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
    * Tumor grade property is not populated <!--SV-585-->
    * Progression_or_recurrence property is not populated <!--SV-584-->
* TARGET projects
    * TARGET CGI BAMs in the Legacy Archive for the following aliquots should not be used because they were not repaired and concatenated into their original composite BAM files by CGHub.
        * TARGET-20-PASJGZ-04A-02D
        * TARGET-30-PAPTLY-01A-01D
        * TARGET-20-PAEIKD-09A-01D
        * TARGET-20-PASMYS-14A-02D
        * TARGET-20-PAMYAS-14A-02D
        * TARGET-10-PAPZST-09A-01D
    * 11 bam files for TARGET-NBL RNA-Seq are not available in the GDC Data portal <!--DAT-1476-->
    * There are 5051 TARGET files for which `experimental_strategy`, `data_format`, `platform`, and `data_subtype` are blank <!--SV-944-->
    * There are two cases with identical submitter_id `TARGET-10-PARUYU` <!--SV-940-->
    * Some TARGET cases are missing `days_to_last_follow_up` <!--SV-934-->
    * Some TARGET cases are missing `age_at_diagnosis` <!--SV-933-->
    * Some TARGET files are not connected to all related aliquots <!--SV-929-->
    * Samples of TARGET sample_type `Recurrent Blood Derived Cancer - Bone Marrow` are mislabeled as `Recurrent Blood Derived Cancer - Peripheral Blood`.  A workaround is to look at the sample barcode, which is -04 for `Recurrent Blood Derived Cancer - Bone Marrow`. (e.g. `TARGET-20-PAMYAS-04A-03R`) <!--SV-918-->
    * The latest TARGET data is not yet available at the GDC.  For the complete and latest data, please see the [TARGET Data Matrix](https://ocg.cancer.gov/programs/target/data-matrix).  Data that is not present or is not the most up to date includes:
        *  All microarray data and metadata
        *  All sequencing analyzed data and metadata
        *  1180 of 12063 sequencing runs of raw data
    * Demographic information for some TARGET patients is incorrect.  The correct information can be found in the associated clinical supplement file.  Impacted patients are TARGET-50-PAJNUS. <!--SV-710-->
    * No data from TARGET-MDLS is available.
* Issues in the Legacy Archive
    * The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
    * SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
    * Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
    * SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
    * TARGET-MDLS cases do not have disease_type or primary_site populated <!--SV-939-->


## Data Release 19.0 <!--REQ-389-->

* __GDC Product__: Data
* __Release Date__: September 17, 2019

### New updates

1.  New projects released:
    *  BEATAML1.0-COHORT - Functional Genomic Landscape of Acute Myeloid Leukemia (phs001657)
        *  Includes WXS and RNA-Seq
2.  New TARGET data released
    *  TARGET-ALL-P1 RNA-Seq
    *  TARGET-ALL-P2 RNA-Seq, WXS, and miRNA-Seq
    *  TARGET-ALL-P3 miRNA-Seq
    *  TARGET-AML WXS, WGS, and miRNA-Seq
    *  TARGET-NBL WXS and RNA-Seq
    *  TARGET-RT WGS and RNA-Seq
    *  TARGET-WT WGS, WXS, and RNA-Seq
3.  Additional CGCI-BLGSP WGS data released
4.  Pindel VCFs released for TARGET-ALL-P2, TARGET-ALL-P3, TARGET-AML, TARGET-NBL, TARGET-WT, MMRF-COMMPASS, HCMI-CMDC, and CPTAC-3
5.  Disease-specific staging properties for many projects were released <!--DAT-2480,  DAT-2503-->


A complete list of files for DR19.0 are listed for the GDC Data Portal and the GDC Legacy Archive are found below:

* [gdc_manifest_20190917_data_release_19.0_active.txt.gz](gdc_manifest_20190917_data_release_19.0_active.txt.gz)
* [gdc_manifest_20190917_data_release_19.0_legacy.txt.gz](gdc_manifest_20190917_data_release_19.0_legacy.txt.gz)


### Bugs Fixed Since Last Release

*  None

### Known Issues and Workarounds

* Some tumor-only annotated VCFs (not raw VCFs) could have a small proportion of variants that appear twice.  Tumor-only annotated VCFs can be identified by searching for workflow "GATK4 MuTect2 Annotation" <!--SV-1425-->
* The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.
* Some miRNA files with QC failed reads were not swapped in DR11.0.  361 aliquots remain to be swapped in a later release <!--DAT-1589-->
* Mutation frequency may be underestimated when using MAF files for genes that overlap other genes.  This is because MAF files only record one gene per variant.
* Most intronic mutations are removed for MAF generation.  However, validated variants may rescue these in some cases.  Therefore intronic mutations in MAF files are not representative of those called by mutation callers.
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
* TCGA Projects
    * Incorrect information about treatment may be included for patients within TCGA-HNSC and TCGA-LGG.  Please refer to the clinical XML for accurate information on treatment <!--DAT-2264, DAT-2265-->
    * 74 Diagnostic TCGA slides are attached to a portion rather than a sample like the rest of the diagnostic slides. The reflects how these original samples were handled. <!--SV-1111-->
    * Two tissue slide images are unavailable for download from GDC Data Portal <!--DAT-1439-->
    * The raw and annotated VarScan VCF files for aliquot `TCGA-VR-A8ET-01A-11D-A403-09` are not available. These VCFs files will be replaced in a later release.<!--TT-602, DAT-1489-->
    * Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
    * Tumor grade property is not populated <!--SV-585-->
    * Progression_or_recurrence property is not populated <!--SV-584-->
* TARGET projects
    * TARGET CGI BAMs in the Legacy Archive for the following aliquots should not be used because they were not repaired and concatenated into their original composite BAM files by CGHub.
        * TARGET-20-PASJGZ-04A-02D
        * TARGET-30-PAPTLY-01A-01D
        * TARGET-20-PAEIKD-09A-01D
        * TARGET-20-PASMYS-14A-02D
        * TARGET-20-PAMYAS-14A-02D
        * TARGET-10-PAPZST-09A-01D
    * 11 bam files for TARGET-NBL RNA-Seq are not available in the GDC Data portal <!--DAT-1476-->
    * There are 5051 TARGET files for which `experimental_strategy`, `data_format`, `platform`, and `data_subtype` are blank <!--SV-944-->
    * There are two cases with identical submitter_id `TARGET-10-PARUYU` <!--SV-940-->
    * Some TARGET cases are missing `days_to_last_follow_up` <!--SV-934-->
    * Some TARGET cases are missing `age_at_diagnosis` <!--SV-933-->
    * Some TARGET files are not connected to all related aliquots <!--SV-929-->
    * Samples of TARGET sample_type `Recurrent Blood Derived Cancer - Bone Marrow` are mislabeled as `Recurrent Blood Derived Cancer - Peripheral Blood`.  A workaround is to look at the sample barcode, which is -04 for `Recurrent Blood Derived Cancer - Bone Marrow`. (e.g. `TARGET-20-PAMYAS-04A-03R`) <!--SV-918-->
    * The latest TARGET data is not yet available at the GDC.  For the complete and latest data, please see the [TARGET Data Matrix](https://ocg.cancer.gov/programs/target/data-matrix).  Data that is not present or is not the most up to date includes:
        *  All microarray data and metadata
        *  All sequencing analyzed data and metadata
        *  1180 of 12063 sequencing runs of raw data
    * Demographic information for some TARGET patients is incorrect.  The correct information can be found in the associated clinical supplement file.  Impacted patients are TARGET-50-PAJNUS. <!--SV-710-->
    * No data from TARGET-MDLS is available.
* Issues in the Legacy Archive
    * The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
    * SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
    * Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
    * SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
    * TARGET-MDLS cases do not have disease_type or primary_site populated <!--SV-939-->


## Data Release 18.0

* __GDC Product__: Data
* __Release Date__: July 8, 2019

### New updates

1.  New Projects released
    *   MMRF-COMMPASS - Multiple Myeloma CoMMpass Study (phs000748)
        *  Includes WGS, WXS, and RNA-Seq
    *   ORGANOID-PANCREATIC - Pancreas Cancer Organoid Profiling (phs001611)
        *  Includes WGS, WXS, and RNA-Seq
    *   TARGET-ALL-P1 - Acute Lymphoblastic Leukemia - Phase I (phs000218)
        *  Includes WGS
    *   TARGET-ALL-P2 - Acute Lymphoblastic Leukemia - Phase II (phs000218)
        *  Includes WGS
    *   CGCI-BLGSP - Burkitt Lymphoma Genome Sequencing Project (phs000235)
        *  Includes WGS and RNA-Seq
2.  New versions of RNA-Seq data for TARGET-ALL-P3     
3.  New RNA-Seq data for TARGET-CCSK
4.  New RNA-Seq data for TARGET-OS


A complete list of files for DR18.0 are listed for the GDC Data Portal and the GDC Legacy Archive are found below:

* [gdc_manifest_20190708_data_release_18.0_active.txt.gz](gdc_manifest_20190708_data_release_18.0_active.txt.gz)
* [gdc_manifest_20190708_data_release_18.0_legacy.txt.gz](gdc_manifest_20190708_data_release_18.0_legacy.txt.gz)


### Bugs Fixed Since Last Release

*  New versions of RNA-Seq data for TARGET-ALL-P3 resolve issue with missing reads from BAM files.

### Known Issues and Workarounds

* Some tumor-only annotated VCFs (not raw VCFs) could have a small proportion of variants that appear twice.  Tumor-only annotated VCFs can be identified by searching for workflow "GATK4 MuTect2 Annotation" <!--SV-1425-->
* The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.
* Some miRNA files with QC failed reads were not swapped in DR11.0.  361 aliquots remain to be swapped in a later release <!--DAT-1589-->
* Mutation frequency may be underestimated when using MAF files for genes that overlap other genes.  This is because MAF files only record one gene per variant.
* Most intronic mutations are removed for MAF generation.  However, validated variants may rescue these in some cases.  Therefore intronic mutations in MAF files are not representative of those called by mutation callers.
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
* TCGA Projects
    * Incorrect information about treatment may be included for patients within TCGA-HNSC and TCGA-LGG.  Please refer to the clinical XML for accurate information on treatment <!--DAT-2264, DAT-2265-->
    * 74 Diagnostic TCGA slides are attached to a portion rather than a sample like the rest of the diagnostic slides. The reflects how these original samples were handled. <!--SV-1111-->
    * Two tissue slide images are unavailable for download from GDC Data Portal <!--DAT-1439-->
    * The raw and annotated VarScan VCF files for aliquot `TCGA-VR-A8ET-01A-11D-A403-09` are not available. These VCFs files will be replaced in a later release.<!--TT-602, DAT-1489-->
    * Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
    * Tumor grade property is not populated <!--SV-585-->
    * Progression_or_recurrence property is not populated <!--SV-584-->
* TARGET projects
    * TARGET CGI BAMs in the Legacy Archive for the following aliquots should not be used because they were not repaired and concatenated into their original composite BAM files by CGHub.
        * TARGET-20-PASJGZ-04A-02D
        * TARGET-30-PAPTLY-01A-01D
        * TARGET-20-PAEIKD-09A-01D
        * TARGET-20-PASMYS-14A-02D
        * TARGET-20-PAMYAS-14A-02D
        * TARGET-10-PAPZST-09A-01D
    * 11 bam files for TARGET-NBL RNA-Seq are not available in the GDC Data portal <!--DAT-1476-->
    * There are 5051 TARGET files for which `experimental_strategy`, `data_format`, `platform`, and `data_subtype` are blank <!--SV-944-->
    * There are two cases with identical submitter_id `TARGET-10-PARUYU` <!--SV-940-->
    * Some TARGET cases are missing `days_to_last_follow_up` <!--SV-934-->
    * Some TARGET cases are missing `age_at_diagnosis` <!--SV-933-->
    * Some TARGET files are not connected to all related aliquots <!--SV-929-->
    * Samples of TARGET sample_type `Recurrent Blood Derived Cancer - Bone Marrow` are mislabeled as `Recurrent Blood Derived Cancer - Peripheral Blood`.  A workaround is to look at the sample barcode, which is -04 for `Recurrent Blood Derived Cancer - Bone Marrow`. (e.g. `TARGET-20-PAMYAS-04A-03R`) <!--SV-918-->
    * The latest TARGET data is not yet available at the GDC.  For the complete and latest data, please see the [TARGET Data Matrix](https://ocg.cancer.gov/programs/target/data-matrix).  Data that is not present or is not the most up to date includes:
        *  All microarray data and metadata
        *  All sequencing analyzed data and metadata
        *  1180 of 12063 sequencing runs of raw data
    * Demographic information for some TARGET patients is incorrect.  The correct information can be found in the associated clinical supplement file.  Impacted patients are TARGET-50-PAJNUS. <!--SV-710-->
    * No data from TARGET-MDLS is available.
* Issues in the Legacy Archive
    * The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
    * SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
    * Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
    * SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
    * TARGET-MDLS cases do not have disease_type or primary_site populated <!--SV-939-->


## Data Release 17.1

* __GDC Product__: Data
* __Release Date__: June 12, 2019

### New updates

1.  Rebuilt indices for NCICCR-DLBCL and CTSP-DLBCL1.  Fewer files viewable in GDC Data Portal or API.

A complete list of files for DR17.1 are listed for the GDC Data Portal and the GDC Legacy Archive are found below:

* [gdc_manifest_20190612_data_release_17.1_active.txt.gz](gdc_manifest_20190612_data_release_17.1_active.txt.gz)
* [gdc_manifest_20190612_data_release_17.1_legacy.txt.gz](gdc_manifest_20190612_data_release_17.1_legacy.txt.gz)

### Bugs Fixed Since Last Release

*  None

### Known Issues and Workarounds


* The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.
* Some miRNA files with QC failed reads were not swapped in DR11.0.  361 aliquots remain to be swapped in a later release <!--DAT-1589-->
* Mutation frequency may be underestimated when using MAF files for genes that overlap other genes.  This is because MAF files only record one gene per variant.
* Most intronic mutations are removed for MAF generation.  However, validated variants may rescue these in some cases.  Therefore intronic mutations in MAF files are not representative of those called by mutation callers.
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
* TCGA Projects
    * Incorrect information about treatment may be included for patients within TCGA-HNSC and TCGA-LGG.  Please refer to the clinical XML for accurate information on treatment <!--DAT-2264, DAT-2265-->
    * 74 Diagnostic TCGA slides are attached to a portion rather than a sample like the rest of the diagnostic slides. The reflects how these original samples were handled. <!--SV-1111-->
    * Two tissue slide images are unavailable for download from GDC Data Portal <!--DAT-1439-->
    * The raw and annotated VarScan VCF files for aliquot `TCGA-VR-A8ET-01A-11D-A403-09` are not available. These VCFs files will be replaced in a later release.<!--TT-602, DAT-1489-->
    * Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
    * Tumor grade property is not populated <!--SV-585-->
    * Progression_or_recurrence property is not populated <!--SV-584-->
* TARGET projects
    * TARGET ALL-P3 RNA-Seq results from DR14 are missing ~18% of reads.  Downsampling appears to be completely random and count files have a very high correlation (>99.99%) with complete data. New versions of these files will be created that include the entire set of reads.
    * TARGET CGI BAMs in the Legacy Archive for the following aliquots should not be used because they were not repaired and concatenated into their original composite BAM files by CGHub.
        * TARGET-20-PASJGZ-04A-02D
        * TARGET-30-PAPTLY-01A-01D
        * TARGET-20-PAEIKD-09A-01D
        * TARGET-20-PASMYS-14A-02D
        * TARGET-20-PAMYAS-14A-02D
        * TARGET-10-PAPZST-09A-01D
    * 11 bam files for TARGET-NBL RNA-Seq are not available in the GDC Data portal <!--DAT-1476-->
    * There are 5051 TARGET files for which `experimental_strategy`, `data_format`, `platform`, and `data_subtype` are blank <!--SV-944-->
    * There are two cases with identical submitter_id `TARGET-10-PARUYU` <!--SV-940-->
    * Some TARGET cases are missing `days_to_last_follow_up` <!--SV-934-->
    * Some TARGET cases are missing `age_at_diagnosis` <!--SV-933-->
    * Some TARGET files are not connected to all related aliquots <!--SV-929-->
    * Samples of TARGET sample_type `Recurrent Blood Derived Cancer - Bone Marrow` are mislabeled as `Recurrent Blood Derived Cancer - Peripheral Blood`.  A workaround is to look at the sample barcode, which is -04 for `Recurrent Blood Derived Cancer - Bone Marrow`. (e.g. `TARGET-20-PAMYAS-04A-03R`) <!--SV-918-->
    * The latest TARGET data is not yet available at the GDC.  For the complete and latest data, please see the [TARGET Data Matrix](https://ocg.cancer.gov/programs/target/data-matrix).  Data that is not present or is not the most up to date includes:
        *  All microarray data and metadata
        *  All sequencing analyzed data and metadata
        *  1180 of 12063 sequencing runs of raw data
    * Demographic information for some TARGET patients is incorrect.  The correct information can be found in the associated clinical supplement file.  Impacted patients are TARGET-50-PAJNUS. <!--SV-710-->
    * No data from TARGET-MDLS is available.
* Issues in the Legacy Archive
    * Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
    * SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
    * Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
    * SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
    * TARGET-MDLS cases do not have disease_type or primary_site populated <!--SV-939-->



## Data Release 17.0 <!--REQ-383-->

* __GDC Product__: Data
* __Release Date__: June 5, 2019

### New updates

1.  New Projects released
    *   HCMI-CMDC - NCI Cancer Model Development for the Human Cancer Model Initiative (HCMI) (phs001486)
    *   BEATAML1.0-CRENOLANIB - Clinical Resistance to Crenolanib in Acute Myeloid Leukemia Due to Diverse Molecular Mechanisms (phs001628)
2.  RNA-Seq data for NCICCR-DLBCL and CTSP-DLBCL1 are released
3.  ATAC-Seq data for TCGA projects are released
4.  CPTAC-3 RNA-Seq data are released
5.  Clinical data updates for TCGA - to see parser code updates review [API v1.20 release notes](https://docs.gdc.cancer.gov/API/Release_Notes/API_Release_Notes/#v1200)
6.  Clinical data updates for other projects to accommodate migration of vital_status, days_to_birth, and days_to_death from the Diagnosis to the Demographic node

A complete list of files for DR17.0 are listed for the GDC Data Portal and the GDC Legacy Archive are found below:

* [gdc_manifest_20190605_data_release_17.0_active.txt.gz](gdc_manifest_20190605_data_release_17.0_active.txt.gz)
* [gdc_manifest_20190605_data_release_17.0_legacy.txt.gz](gdc_manifest_20190605_data_release_17.0_legacy.txt.gz).

### Bugs Fixed Since Last Release

*  None

### Known Issues and Workarounds


* The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.
* Some miRNA files with QC failed reads were not swapped in DR11.0.  361 aliquots remain to be swapped in a later release <!--DAT-1589-->
* Mutation frequency may be underestimated when using MAF files for genes that overlap other genes.  This is because MAF files only record one gene per variant.
* Most intronic mutations are removed for MAF generation.  However, validated variants may rescue these in some cases.  Therefore intronic mutations in MAF files are not representative of those called by mutation callers.
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
* TCGA Projects
    * Incorrect information about treatment may be included for patients within TCGA-HNSC and TCGA-LGG.  Please refer to the clinical XML for accurate information on treatment <!--DAT-2264, DAT-2265-->
    * 74 Diagnostic TCGA slides are attached to a portion rather than a sample like the rest of the diagnostic slides. The reflects how these original samples were handled. <!--SV-1111-->
    * Two tissue slide images are unavailable for download from GDC Data Portal <!--DAT-1439-->
    * The raw and annotated VarScan VCF files for aliquot `TCGA-VR-A8ET-01A-11D-A403-09` are not available. These VCFs files will be replaced in a later release.<!--TT-602, DAT-1489-->
    * Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
    * Tumor grade property is not populated <!--SV-585-->
    * Progression_or_recurrence property is not populated <!--SV-584-->
* TARGET projects
    * TARGET ALL-P3 RNA-Seq results from DR14 are missing ~18% of reads.  Downsampling appears to be completely random and count files have a very high correlation (>99.99%) with complete data. New versions of these files will be created that include the entire set of reads.
    * TARGET CGI BAMs in the Legacy Archive for the following aliquots should not be used because they were not repaired and concatenated into their original composite BAM files by CGHub.
        * TARGET-20-PASJGZ-04A-02D
        * TARGET-30-PAPTLY-01A-01D
        * TARGET-20-PAEIKD-09A-01D
        * TARGET-20-PASMYS-14A-02D
        * TARGET-20-PAMYAS-14A-02D
        * TARGET-10-PAPZST-09A-01D
    * 11 bam files for TARGET-NBL RNA-Seq are not available in the GDC Data portal <!--DAT-1476-->
    * There are 5051 TARGET files for which `experimental_strategy`, `data_format`, `platform`, and `data_subtype` are blank <!--SV-944-->
    * There are two cases with identical submitter_id `TARGET-10-PARUYU` <!--SV-940-->
    * Some TARGET cases are missing `days_to_last_follow_up` <!--SV-934-->
    * Some TARGET cases are missing `age_at_diagnosis` <!--SV-933-->
    * Some TARGET files are not connected to all related aliquots <!--SV-929-->
    * Samples of TARGET sample_type `Recurrent Blood Derived Cancer - Bone Marrow` are mislabeled as `Recurrent Blood Derived Cancer - Peripheral Blood`.  A workaround is to look at the sample barcode, which is -04 for `Recurrent Blood Derived Cancer - Bone Marrow`. (e.g. `TARGET-20-PAMYAS-04A-03R`) <!--SV-918-->
    * The latest TARGET data is not yet available at the GDC.  For the complete and latest data, please see the [TARGET Data Matrix](https://ocg.cancer.gov/programs/target/data-matrix).  Data that is not present or is not the most up to date includes:
        *  All microarray data and metadata
        *  All sequencing analyzed data and metadata
        *  1180 of 12063 sequencing runs of raw data
    * Demographic information for some TARGET patients is incorrect.  The correct information can be found in the associated clinical supplement file.  Impacted patients are TARGET-50-PAJNUS. <!--SV-710-->
    * No data from TARGET-MDLS is available.
* Issues in the Legacy Archive
    * Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
    * SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
    * Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
    * SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
    * TARGET-MDLS cases do not have disease_type or primary_site populated <!--SV-939-->





## Data Release 16.0 <!--REQ-385-->

* __GDC Product__: Data
* __Release Date__: March 26, 2019

### New updates

1.  The CPTAC-3 project (phs001287) is released with WXS and WGS data.  RNA-Seq will be released at a later date.  Additional project details can be found at on the [CPTAC Data Source page](https://gdc.cancer.gov/about-gdc/contributed-genomic-data-cancer-research/clinical-proteomic-tumor-analysis-consortium-cptac).
2.  TARGET-ALL-P3 (phs000218) WGS BAM files are released.
3.  VAREPOP-APOLLO (phs001374) VCF files are released.

A complete list of files for DR16.0 are listed for the GDC Data Portal and the GDC Legacy Archive are found below:

* [gdc_manifest_20190326_data_release_16.0_active.txt.gz](gdc_manifest_20190326_data_release_16.0_active.txt.gz)
* [gdc_manifest_20190326_data_release_16.0_legacy.txt.gz](gdc_manifest_20190326_data_release_16.0_legacy.txt.gz).

### Bugs Fixed Since Last Release

*  None

### Known Issues and Workarounds

* TARGET ALL-P3 RNA-Seq results from DR14 are missing ~18% of reads.  Downsampling appears to be completely random and count files have a very high correlation (>99.99%) with complete data. New versions of these files will be created that include the entire set of reads.
* The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.
* TARGET CGI BAMs in the Legacy Archive for the following aliquots should not be used because they were not repaired and concatenated into their original composite BAM files by CGHub.
    * TARGET-20-PASJGZ-04A-02D
    * TARGET-30-PAPTLY-01A-01D
    * TARGET-20-PAEIKD-09A-01D
    * TARGET-20-PASMYS-14A-02D
    * TARGET-20-PAMYAS-14A-02D
    * TARGET-10-PAPZST-09A-01D
* Some miRNA files with QC failed reads were not swapped in DR11.0.  361 aliquots remain to be swapped in a later release <!--DAT-1589-->
* 74 Diagnostic TCGA slides are attached to a portion rather than a sample like the rest of the diagnostic slides. The reflects how these original samples were handled. <!--SV-1111-->
* 11 bam files for TARGET-NBL RNA-Seq are not available in the GDC Data portal <!--DAT-1476-->
* Two tissue slide images are unavailable for download from GDC Data Portal <!--DAT-1439-->
* The raw and annotated VarScan VCF files for aliquot `TCGA-VR-A8ET-01A-11D-A403-09` are not available. These VCFs files will be replaced in a later release.
<!--TT-602, DAT-1489-->
* There are 5051 TARGET files for which `experimental_strategy`, `data_format`, `platform`, and `data_subtype` are blank <!--SV-944-->
* There are two cases with identical submitter_id `TARGET-10-PARUYU` <!--SV-940-->
* TARGET-MDLS cases do not have disease_type or primary_site populated <!--SV-939-->
* Some TARGET cases are missing `days_to_last_follow_up` <!--SV-934-->
* Some TARGET cases are missing `age_at_diagnosis` <!--SV-933-->
* Some TARGET files are not connected to all related aliquots <!--SV-929-->
* Samples of TARGET sample_type `Recurrent Blood Derived Cancer - Bone Marrow` are mislabeled as `Recurrent Blood Derived Cancer - Peripheral Blood`.  A workaround is to look at the sample barcode, which is -04 for `Recurrent Blood Derived Cancer - Bone Marrow`. (e.g. `TARGET-20-PAMYAS-04A-03R`) <!--SV-918-->
* Mutation frequency may be underestimated when using MAF files for genes that overlap other genes.  This is because MAF files only record one gene per variant.
* Most intronic mutations are removed for MAF generation.  However, validated variants may rescue these in some cases.  Therefore intronic mutations in MAF files are not representative of those called by mutation callers.
* The latest TARGET data is not yet available at the GDC.  For the complete and latest data, please see the [TARGET Data Matrix](https://ocg.cancer.gov/programs/target/data-matrix).  Data that is not present or is not the most up to date includes:
    *  All microarray data and metadata
    *  All sequencing analyzed data and metadata
    *  1180 of 12063 sequencing runs of raw data
* Demographic information for some TARGET patients is incorrect.  The correct information can be found in the associated clinical supplement file.  Impacted patients are TARGET-50-PAJNUS. <!--SV-710-->
* Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* No data from TARGET-MDLS is available.
* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
* SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
* Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
* SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
* Tumor grade property is not populated <!--SV-585-->
* Progression_or_recurrence property is not populated <!--SV-584-->




## Data Release 15.0 <!--REQ-381-->

* __GDC Product__: Data
* __Release Date__: February 20, 2019

### New updates

1.  TARGET-ALL-P3 is now available and includes RNA-Seq and WXS data.
2.  New RNA-Seq workflow is now being utilized for new projects.  More details can be found in the [RNA-Seq pipeline documentation](../../Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/#rna-seq-alignment-workflow).
3.  New tumor only variant calling pipeline is now being utilized for new projects.  More details can be found in the [Tumor only pipeline documentation](../../Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/#tumor-only-variant-calling-workflow).


A complete list of files for DR15.0 are listed for the GDC Data Portal and the GDC Legacy Archive are found below:

* [gdc_manifest_20190220_data_release_15.0_active.txt.gz](gdc_manifest_20190220_data_release_15.0_active.txt.gz)
* [gdc_manifest_20190220_data_release_15.0_legacy.txt.gz](gdc_manifest_20190220_data_release_15.0_legacy.txt.gz).

### Bugs Fixed Since Last Release

*  None

### Known Issues and Workarounds

* The read alignment end coordinates in the x.isoform.quantification.txt files produced by the miRNA pipeline are exclusive (i.e. offset by 1) for all TCGA miRNA legacy (GRCh37/hg19) and current harmonized (GRCh38/hg38) miRNA data.  This error has no impact on miRNA alignment or quantification - only the coordinates reported in the quantification file.
* TARGET CGI BAMs in the Legacy Archive for the following aliquots should not be used because they were not repaired and concatenated into their original composite BAM files by CGHub.
    * TARGET-20-PASJGZ-04A-02D
    * TARGET-30-PAPTLY-01A-01D
    * TARGET-20-PAEIKD-09A-01D
    * TARGET-20-PASMYS-14A-02D
    * TARGET-20-PAMYAS-14A-02D
    * TARGET-10-PAPZST-09A-01D
* Some miRNA files with QC failed reads were not swapped in DR11.0.  361 aliquots remain to be swapped in a later release <!--DAT-1589-->
* 74 Diagnostic TCGA slides are attached to a portion rather than a sample like the rest of the diagnostic slides. The reflects how these original samples were handled. <!--SV-1111-->
* 11 bam files for TARGET-NBL RNA-Seq are not available in the GDC Data portal <!--DAT-1476-->
* Two tissue slide images are unavailable for download from GDC Data Portal <!--DAT-1439-->
* The raw and annotated VarScan VCF files for aliquot `TCGA-VR-A8ET-01A-11D-A403-09` are not available. These VCFs files will be replaced in a later release.
<!--TT-602, DAT-1489-->
* There are 5051 TARGET files for which `experimental_strategy`, `data_format`, `platform`, and `data_subtype` are blank <!--SV-944-->
* There are two cases with identical submitter_id `TARGET-10-PARUYU` <!--SV-940-->
* TARGET-MDLS cases do not have disease_type or primary_site populated <!--SV-939-->
* Some TARGET cases are missing `days_to_last_follow_up` <!--SV-934-->
* Some TARGET cases are missing `age_at_diagnosis` <!--SV-933-->
* Some TARGET files are not connected to all related aliquots <!--SV-929-->
* Samples of TARGET sample_type `Recurrent Blood Derived Cancer - Bone Marrow` are mislabeled as `Recurrent Blood Derived Cancer - Peripheral Blood`.  A workaround is to look at the sample barcode, which is -04 for `Recurrent Blood Derived Cancer - Bone Marrow`. (e.g. `TARGET-20-PAMYAS-04A-03R`) <!--SV-918-->
* Mutation frequency may be underestimated when using MAF files for genes that overlap other genes.  This is because MAF files only record one gene per variant.
* Most intronic mutations are removed for MAF generation.  However, validated variants may rescue these in some cases.  Therefore intronic mutations in MAF files are not representative of those called by mutation callers.
* The latest TARGET data is not yet available at the GDC.  For the complete and latest data, please see the [TARGET Data Matrix](https://ocg.cancer.gov/programs/target/data-matrix).  Data that is not present or is not the most up to date includes:
    *  All microarray data and metadata
    *  All sequencing analyzed data and metadata
    *  1180 of 12063 sequencing runs of raw data
* Demographic information for some TARGET patients is incorrect.  The correct information can be found in the associated clinical supplement file.  Impacted patients are TARGET-50-PAJNUS. <!--SV-710-->
* Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* No data from TARGET-MDLS is available.
* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
* SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
* Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
* SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
* Tumor grade property is not populated <!--SV-585-->
* Progression_or_recurrence property is not populated <!--SV-584-->


## Data Release 14.0 <!--REQ-335-->

* __GDC Product__: Data
* __Release Date__: December 18, 2018

### New updates

1.  Copy Number Variation (CNV) data derived from GISTIC2 results are now available for download for TCGA projects <!--DAT-1681-->
2.  New miRNA data available for 181 aliquots for TARGET and TCGA <!--DAT-1823-->
3.  Released two SNP6 files (6cd4ef5e-324a-4ace-8779-7a33bd559c83, dfa89ee9-6ee5-460b-bd58-b5ca0e9cb7ac) <!--DAT-1805-->
4.  New versions of TCGA biospecimen supplements are available <!--DAT-1776-->
5.  Updated primary site for `TCGA-AG-3881` to `Unknown` <!--DAT-1853-->
6.  8 New Harmonized WGS BAM files for TARGET-WT, TARGET-NBL, TARGET-AML added to the portal <!--DAT-1839-->

A complete list of files for DR14.0 are listed for the GDC Data Portal and the GDC Legacy Archive are found below:

* [gdc_manifest_20181218_data_release_14.0_active.txt.gz](gdc_manifest_20181218_data_release_14.0_active.txt.gz)
* [gdc_manifest_20181218_data_release_14.0_legacy.txt.gz](gdc_manifest_20181218_data_release_14.0_legacy.txt.gz).


### Bugs Fixed Since Last Release

*  FM-AD clinial and biospecimen supplements are now correctly labeled as TSV rather than XLSX <!--DAT-1681-->

### Known Issues and Workarounds

* TARGET CGI BAMs in the Legacy Archive for the following aliquots should not be used because they were not repaired and concatenated into their original composite BAM files by CGHub.
    * TARGET-20-PASJGZ-04A-02D
    * TARGET-30-PAPTLY-01A-01D
    * TARGET-20-PAEIKD-09A-01D
    * TARGET-20-PASMYS-14A-02D
    * TARGET-20-PAMYAS-14A-02D
    * TARGET-10-PAPZST-09A-01D
* Some miRNA files with QC failed reads were not swapped in DR11.0.  361 aliquots remain to be swapped in a later release <!--DAT-1589-->
* 74 Diagnostic TCGA slides are attached to a portion rather than a sample like the rest of the diagnostic slides. The reflects how these original samples were handled. <!--SV-1111-->
* 11 bam files for TARGET-NBL RNA-Seq are not available in the GDC Data portal <!--DAT-1476-->
* Two tissue slide images are unavailable for download from GDC Data Portal <!--DAT-1439-->
* The raw and annotated VarScan VCF files for aliquot `TCGA-VR-A8ET-01A-11D-A403-09` are not available. These VCFs files will be replaced in a later release.
<!--TT-602, DAT-1489-->
* There are 5051 TARGET files for which `experimental_strategy`, `data_format`, `platform`, and `data_subtype` are blank <!--SV-944-->
* There are two cases with identical submitter_id `TARGET-10-PARUYU` <!--SV-940-->
* TARGET-MDLS cases do not have disease_type or primary_site populated <!--SV-939-->
* Some TARGET cases are missing `days_to_last_follow_up` <!--SV-934-->
* Some TARGET cases are missing `age_at_diagnosis` <!--SV-933-->
* Some TARGET files are not connected to all related aliquots <!--SV-929-->
* Samples of TARGET sample_type `Recurrent Blood Derived Cancer - Bone Marrow` are mislabeled as `Recurrent Blood Derived Cancer - Peripheral Blood`.  A workaround is to look at the sample barcode, which is -04 for `Recurrent Blood Derived Cancer - Bone Marrow`. (e.g. `TARGET-20-PAMYAS-04A-03R`) <!--SV-918-->
* Mutation frequency may be underestimated when using MAF files for genes that overlap other genes.  This is because MAF files only record one gene per variant.
* Most intronic mutations are removed for MAF generation.  However, validated variants may rescue these in some cases.  Therefore intronic mutations in MAF files are not representative of those called by mutation callers.
* The latest TARGET data is not yet available at the GDC.  For the complete and latest data, please see the [TARGET Data Matrix](https://ocg.cancer.gov/programs/target/data-matrix).  Data that is not present or is not the most up to date includes:
    *  All microarray data and metadata
    *  All sequencing analyzed data and metadata
    *  1180 of 12063 sequencing runs of raw data
* Demographic information for some TARGET patients is incorrect.  The correct information can be found in the associated clinical supplement file.  Impacted patients are TARGET-50-PAJNUS. <!--SV-710-->
* Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* No data from TARGET-MDLS is available.
* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
* SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
* Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
* SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
* Tumor grade property is not populated <!--SV-585-->
* Progression_or_recurrence property is not populated <!--SV-584-->



## Data Release 13.0

* __GDC Product__: Data
* __Release Date__: September 27, 2018

### New updates

1. Three new projects are released to the GDC (VAREPOP-APOLLO (phs001374), CTSP-DLBCL1 (phs001184), NCICCR-DLBCL (phs001444) <!--SPT-136,SPT-145,SPT-84-->
2. TARGET WGS alignments are released. VCFs will be provided in a later release <!--DAT-1699-->
3. Clinical data was harmonized with ICD-O-3 terminology for TCGA properties case.primary_site, case.disease_type, diagnosis.primary_diagnosis, diagnosis.site_of_resection_or_biopsy, diagnosis.tissue_or_organ_of_origin <!--DAT-1512-->
4. Redaction annotations applied to 11 aliquots in TCGA-DLBC <!--DAT-1169, Data-871-->
5. Redaction annotations applied to incorrectly trimmed miRNA file in the Legacy Achive <!--SV-653-->

A complete list of files for DR13.0 are listed for the GDC Data Portal and the GDC Legacy Archive are found below:

* [gdc_manifest_20180927_data_release_13.0_active.txt.gz](gdc_manifest_20180927_data_release_13.0_active.txt.gz)
* [gdc_manifest_20180927_data_release_13.0_legacy.txt.gz](gdc_manifest_20180927_data_release_13.0_legacy.txt.gz).


### Bugs Fixed Since Last Release

* 253 files Copy Number Segment and Masked Copy Number Segment files were released.  These were skipped in DR 12.0 <!--DAT-1404-->
* 36 Diagnostic TCGA slides were released.  They were skipped in DR 12.0  <!--SV-1109-->

### Known Issues and Workarounds

* 506 Copy Number Segment and 36 Slide Image files are designated as controlled-access on the GDC Data Portal.  These files are actually open-access and will be downloadable without a token using [this manifest](gdc_manifest_20181003_data_release_13.0_cnv_slides.txt).  <!--DAT-1804-->
* 2 Copy Number Segment files from TCGA-TGCT do not appear on the GDC Portal. They can be downloaded using the Data Transfer Tool using the following UUIDs. <!--DAT-1805-->
    * 6cd4ef5e-324a-4ace-8779-7a33bd559c83 - RAMPS_p_TCGA_Batch_430_NSP_GenomeWideSNP_6_E07_1538238.nocnv_grch38.seg.v2.txt
    * dfa89ee9-6ee5-460b-bd58-b5ca0e9cb7ac - RAMPS_p_TCGA_Batch_430_NSP_GenomeWideSNP_6_E07_1538238.grch38.seg.v2.txt
* TARGET CGI BAMs in the Legacy Archive for the following aliquots should not be used because they were not repaired and concatenated into their original composite BAM files by CGHub.
    * TARGET-20-PASJGZ-04A-02D
    * TARGET-30-PAPTLY-01A-01D
    * TARGET-20-PAEIKD-09A-01D
    * TARGET-20-PASMYS-14A-02D
    * TARGET-20-PAMYAS-14A-02D
    * TARGET-10-PAPZST-09A-01D
* Some miRNA files with QC failed reads were not swapped in DR11.0.  361 aliquots remain to be swapped in a later release <!--DAT-1589-->
* 74 Diagnostic TCGA slides are attached to a portion rather than a sample like the rest of the diagnostic slides. The reflects how these original samples were handled. <!--SV-1111-->
* 11 bam files for TARGET-NBL RNA-Seq are not available in the GDC Data portal <!--DAT-1476-->
* Two tissue slide images are unavailable for download from GDC Data Portal <!--DAT-1439-->
* The raw and annotated VarScan VCF files for aliquot `TCGA-VR-A8ET-01A-11D-A403-09` are not available. These VCFs files will be replaced in a later release.
<!--TT-602, DAT-1489-->
* There are 5051 TARGET files for which `experimental_strategy`, `data_format`, `platform`, and `data_subtype` are blank <!--SV-944-->
* There are two cases with identical submitter_id `TARGET-10-PARUYU` <!--SV-940-->
* TARGET-MDLS cases do not have disease_type or primary_site populated <!--SV-939-->
* Some TARGET cases are missing `days_to_last_follow_up` <!--SV-934-->
* Some TARGET cases are missing `age_at_diagnosis` <!--SV-933-->
* Some TARGET files are not connected to all related aliquots <!--SV-929-->
* Samples of TARGET sample_type `Recurrent Blood Derived Cancer - Bone Marrow` are mislabeled as `Recurrent Blood Derived Cancer - Peripheral Blood`.  A workaround is to look at the sample barcode, which is -04 for `Recurrent Blood Derived Cancer - Bone Marrow`. (e.g. `TARGET-20-PAMYAS-04A-03R`) <!--SV-918-->
* FM-AD clinical and biospecimen supplement files have incorrect data format.  They are listed as XLSX, but are in fact TSV files. <!--DAT-1123-->
* Mutation frequency may be underestimated when using MAF files for genes that overlap other genes.  This is because MAF files only record one gene per variant.
* Most intronic mutations are removed for MAF generation.  However, validated variants may rescue these in some cases.  Therefore intronic mutations in MAF files are not representative of those called by mutation callers.
* The latest TARGET data is not yet available at the GDC.  For the complete and latest data, please see the [TARGET Data Matrix](https://ocg.cancer.gov/programs/target/data-matrix).  Data that is not present or is not the most up to date includes:
    *  All microarray data and metadata
    *  All sequencing analyzed data and metadata
    *  1180 of 12063 sequencing runs of raw data
* Demographic information for some TARGET patients is incorrect.  The correct information can be found in the associated clinical supplement file.  Impacted patients are TARGET-50-PAJNUS. <!--SV-710-->
* Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* No data from TARGET-MDLS is available.
* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
* SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
* Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
* SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
* Tumor grade property is not populated <!--SV-585-->
* Progression_or_recurrence property is not populated <!--SV-584-->




## Data Release 12.0

* __GDC Product__: Data
* __Release Date__: June 13, 2018

### New updates

1. Updated clinical and biospecimen XML files for TCGA cases are available in the GDC Data Portal. Equivalent Legacy Archive files may no longer be up to date. <!--TT-328-->
2. All biospecimen and clinical supplement files for TCGA projects formerly only found in the Legacy Archive have been updated and transferred to the GDC Data Portal. Equivalent Legacy Archive files and metadata retrieved from the API may no longer be up to date. <!--TT-327-->
3. Diagnostic slides from TCGA are now available in the GDC Data Portal and Slide Image Viewer.  They were formerly only available in the Legacy Archive. <!--DAT-1316-->
4. Updated Copy Number Segment and Masked Copy Number Segment files are now available.  These were generated using an improved mapping of hg38 coordinates for the Affymetrix SNP6.0 probe set. <!--DAT-1303-->
5.  VCF files containing SNVs produced from TARGET WGS CGI data are available.  The variant calls were initially produced by CGI and lifted over to hg38. <!--DAT-1281-->

Updated files for this release are listed [here](DR12.0_files_swap.txt.gz).
A complete list of files for DR12.0 are listed for the GDC Data Portal [here](gdc_manifest_20180613_data_release_12.0_active.txt.gz) and the GDC Legacy Archive [here](gdc_manifest_20180613_data_release_12.0_legacy.txt.gz).


### Bugs Fixed Since Last Release

* TARGET NBL RNA-Seq data is now associated with the correct aliquot. <!--SV-1097-->

### Known Issues and Workarounds

* Some Copy Number Segment and Masked Copy Number Segment were not replaced in DR 12.0.  253 files remain to be swapped in a later release <!--DAT-1404-->
* Some miRNA files with QC failed reads were not swapped in DR11.0.  361 aliquots remain to be swapped in a later release <!--DAT-1589-->
* 74 Diagnostic TCGA slides are attached to a portion rather than a sample like the rest of the diagnostic slides. The reflects how these original samples were handled. <!--SV-1111-->
* 36 Diagnostic TCGA slides are not yet available in the active GDC Portal. They are still available in the GDC Legacy Archive.  <!--SV-1109-->
* 11 bam files for TARGET-NBL RNA-Seq are not available in the GDC Data portal <!--DAT-1476-->
* Two tissue slide images are unavailable for download from GDC Data Portal <!--DAT-1439-->
* The raw and annotated VarScan VCF files for aliquot `TCGA-VR-A8ET-01A-11D-A403-09` are not available. These VCFs files will be replaced in a later release.
<!--TT-602, DAT-1489-->
* There are 5051 TARGET files for which `experimental_strategy`, `data_format`, `platform`, and `data_subtype` are blank <!--SV-944-->
* There are two cases with identical submitter_id `TARGET-10-PARUYU` <!--SV-940-->
* TARGET-MDLS cases do not have disease_type or primary_site populated <!--SV-939-->
* Some TARGET cases are missing `days_to_last_follow_up` <!--SV-934-->
* Some TARGET cases are missing `age_at_diagnosis` <!--SV-933-->
* Some TARGET files are not connected to all related aliquots <!--SV-929-->
* Samples of TARGET sample_type `Recurrent Blood Derived Cancer - Bone Marrow` are mislabeled as `Recurrent Blood Derived Cancer - Peripheral Blood`.  A workaround is to look at the sample barcode, which is -04 for `Recurrent Blood Derived Cancer - Bone Marrow`. (e.g. `TARGET-20-PAMYAS-04A-03R`) <!--SV-918-->
* FM-AD clinical and biospecimen supplement files have incorrect data format.  They are listed as XLSX, but are in fact TSV files. <!--DAT-1123-->
* Mutation frequency may be underestimated when using MAF files for genes that overlap other genes.  This is because MAF files only record one gene per variant.
* Most intronic mutations are removed for MAF generation.  However, validated variants may rescue these in some cases.  Therefore intronic mutations in MAF files are not representative of those called by mutation callers.
* The latest TARGET data is not yet available at the GDC.  For the complete and latest data, please see the [TARGET Data Matrix](https://ocg.cancer.gov/programs/target/data-matrix).  Data that is not present or is not the most up to date includes:
    *  All microarray data and metadata
    *  All sequencing analyzed data and metadata
    *  1180 of 12063 sequencing runs of raw data
* Demographic information for some TARGET patients is incorrect.  The correct information can be found in the associated clinical supplement file.  Impacted patients are TARGET-50-PAJNUS. <!--SV-710-->
* There are 11 cases in project TCGA-DLBC that are known to have incorrect WXS data in the GDC Data Portal.  Impacted cases are listed below.  This affects the BAMs and VCFs associated with these cases in the GDC Data Portal.  Corrected BAMs can be found in the GDC Legacy Archive.  Variants from affected aliquots appear in the protected MAFs with GDC_FILTER=ContEst to indicate a sample contamination problem, but are removed during the generation of the Somatic MAF file.  In a later release we will supply corrected BAM, VCF, and MAF files for these cases.  In the mean time, we advise you not to use any of the WXS files associated with these cases in the GDC Data Portal.  A list of these files can be found here. [Download list of affected files](DLBC_Affected_Files.txt). <!-- Data-871-->
    * TCGA-FF-8062
    * TCGA-FM-8000
    * TCGA-G8-6324
    * TCGA-G8-6325
    * TCGA-G8-6326
    * TCGA-G8-6906
    * TCGA-G8-6907
    * TCGA-G8-6909
    * TCGA-G8-6914
    * TCGA-GR-7351
    * TCGA-GR-7353
* Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* No data from TARGET-MDLS is available.
* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
* SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
* Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
* SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
* Tumor grade property is not populated <!--SV-585-->
* Progression_or_recurrence property is not populated <!--SV-584-->




## Data Release 11.0

* __GDC Product__: Data
* __Release Date__: May 21, 2018

### New updates

1. Updated miRNA files to remove QCFail reads.  This included all BAM and downstream count files. <!--DAT-1246-->
2. TCGA Tissue slide images now available in GDC Data Portal.  Previously these were found only in the Legacy Archive <!--DAT-1251-->

Updated files for this release are listed [here](DR11.0_files_swap.txt.gz).
A complete list of files for DR11.0 are listed for the GDC Data Portal [here](gdc_manifest_20180521_data_release_11.0_active.txt.gz) and the GDC Legacy Archive [here](gdc_manifest_20180521_data_release_11.0_legacy.txt.gz).


### Bugs Fixed Since Last Release

* N/A

### Known Issues and Workarounds

* Two tissue slide images are unavailable for download from GDC Data Portal <!--DAT-1439-->
* RNA-Seq files for TARGET-NBL are attached to the incorrect aliquot.  The BAM files contain the correct information in their header but the connection in the GDC to read groups and aliquots is incorrect.  The linked file below contains a mapping between aliquots where file are currently associated and the aliquot where they should instead be associated [(mapping file)](correct_aliquot_mappings.tsv). <!--SV-1097-->
* The raw and annotated VarScan VCF files for aliquot `TCGA-VR-A8ET-01A-11D-A403-09` were not replaced in DR10.0 and thus do not contain indels.  However, the indels from this aliquot can be found in the MAF files and are displayed in the Exploration section in the Data Portal.  These VCFs files will be replaced in a later release.
<!--SV-950-->
* There are 5051 TARGET files for which `experimental_strategy`, `data_format`, `platform`, and `data_subtype` are blank <!--SV-944-->
* There are two cases with identical submitter_id `TARGET-10-PARUYU` <!--SV-940-->
* TARGET-MDLS cases do not have disease_type or primary_site populated <!--SV-939-->
* Some TARGET cases are missing `days_to_last_follow_up` <!--SV-934-->
* Some TARGET cases are missing `age_at_diagnosis` <!--SV-933-->
* Some TARGET files are not connected to all related aliquots <!--SV-929-->
* miRNA alignments include QC failed reads.  
* Samples of TARGET sample_type `Recurrent Blood Derived Cancer - Bone Marrow` are mislabeled as `Recurrent Blood Derived Cancer - Peripheral Blood`.  A workaround is to look at the sample barcode, which is -04 for `Recurrent Blood Derived Cancer - Bone Marrow`. (e.g. `TARGET-20-PAMYAS-04A-03R`) <!--SV-918-->
* FM-AD clinical and biospecimen supplement files have incorrect data format.  They are listed as XLSX, but are in fact TSV files. <!--DAT-1123-->
* Mutation frequency may be underestimated when using MAF files for genes that overlap other genes.  This is because MAF files only record one gene per variant.
* Most intronic mutations are removed for MAF generation.  However, validated variants may rescue these in some cases.  Therefore intronic mutations in MAF files are not representative of those called by mutation callers.
* The latest TARGET data is not yet available at the GDC.  For the complete and latest data, please see the [TARGET Data Matrix](https://ocg.cancer.gov/programs/target/data-matrix).  Data that is not present or is not the most up to date includes:
    *  All microarray data and metadata
    *  All sequencing analyzed data and metadata
    *  1180 of 12063 sequencing runs of raw data
* Demographic information for some TARGET patients is incorrect.  The correct information can be found in the associated clinical supplement file.  Impacted patients are TARGET-50-PAJNUS. <!--SV-710-->
* There are 11 cases in project TCGA-DLBC that are known to have incorrect WXS data in the GDC Data Portal.  Impacted cases are listed below.  This affects the BAMs and VCFs associated with these cases in the GDC Data Portal.  Corrected BAMs can be found in the GDC Legacy Archive.  Variants from affected aliquots appear in the protected MAFs with GDC_FILTER=ContEst to indicate a sample contamination problem, but are removed during the generation of the Somatic MAF file.  In a later release we will supply corrected BAM, VCF, and MAF files for these cases.  In the mean time, we advise you not to use any of the WXS files associated with these cases in the GDC Data Portal.  A list of these files can be found here. [Download list of affected files](DLBC_Affected_Files.txt). <!-- Data-871-->
    * TCGA-FF-8062
    * TCGA-FM-8000
    * TCGA-G8-6324
    * TCGA-G8-6325
    * TCGA-G8-6326
    * TCGA-G8-6906
    * TCGA-G8-6907
    * TCGA-G8-6909
    * TCGA-G8-6914
    * TCGA-GR-7351
    * TCGA-GR-7353
* Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* No data from TARGET-MDLS is available.
* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
* SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
* Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
* SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
* Tumor grade property is not populated <!--SV-585-->
* Progression_or_recurrence property is not populated <!--SV-584-->



## Data Release 10.1

* __GDC Product__: Data
* __Release Date__: February 15, 2018

### New updates

1. Updated FM-AD clinical data to conform with Data Dictionary release v1.11

### Bugs Fixed Since Last Release

None

### Known Issues and Workarounds

* RNA-Seq files for TARGET-NBL are attached to the incorrect aliquot.  The BAM files contain the correct information in their header but the connection in the GDC to read groups and aliquots is incorrect.  The linked file below contains a mapping between aliquots where file are currently associated and the aliquot where they should instead be associated [(mapping file)](correct_aliquot_mappings.tsv). <!--SV-1097-->
* The raw and annotated VarScan VCF files for aliquot `TCGA-VR-A8ET-01A-11D-A403-09` were not replaced in DR10.0 and thus do not contain indels.  However, the indels from this aliquot can be found in the MAF files and are displayed in the Exploration section in the Data Portal.  These VCFs files will be replaced in a later release.
<!--SV-950-->
* There are 5051 TARGET files for which `experimental_strategy`, `data_format`, `platform`, and `data_subtype` are blank <!--SV-944-->
* There are two cases with identical submitter_id `TARGET-10-PARUYU` <!--SV-940-->
* TARGET-MDLS cases do not have disease_type or primary_site populated <!--SV-939-->
* Some TARGET cases are missing `days_to_last_follow_up` <!--SV-934-->
* Some TARGET cases are missing `age_at_diagnosis` <!--SV-933-->
* Some TARGET files are not connected to all related aliquots <!--SV-929-->
* miRNA alignments include QC failed reads.  
* Samples of TARGET sample_type `Recurrent Blood Derived Cancer - Bone Marrow` are mislabeled as `Recurrent Blood Derived Cancer - Peripheral Blood`.  A workaround is to look at the sample barcode, which is -04 for `Recurrent Blood Derived Cancer - Bone Marrow`. (e.g. `TARGET-20-PAMYAS-04A-03R`) <!--SV-918-->
* FM-AD clinical and biospecimen supplement files have incorrect data format.  They are listed as XLSX, but are in fact TSV files. <!--DAT-1123-->
* Mutation frequency may be underestimated when using MAF files for genes that overlap other genes.  This is because MAF files only record one gene per variant.
* Most intronic mutations are removed for MAF generation.  However, validated variants may rescue these in some cases.  Therefore intronic mutations in MAF files are not representative of those called by mutation callers.
* The latest TARGET data is not yet available at the GDC.  For the complete and latest data, please see the [TARGET Data Matrix](https://ocg.cancer.gov/programs/target/data-matrix).  Data that is not present or is not the most up to date includes:
    *  All microarray data and metadata
    *  All sequencing analyzed data and metadata
    *  1180 of 12063 sequencing runs of raw data
* Demographic information for some TARGET patients is incorrect.  The correct information can be found in the associated clinical supplement file.  Impacted patients are TARGET-50-PAJNUS. <!--SV-710-->
* There are 11 cases in project TCGA-DLBC that are known to have incorrect WXS data in the GDC Data Portal.  Impacted cases are listed below.  This affects the BAMs and VCFs associated with these cases in the GDC Data Portal.  Corrected BAMs can be found in the GDC Legacy Archive.  Variants from affected aliquots appear in the protected MAFs with GDC_FILTER=ContEst to indicate a sample contamination problem, but are removed during the generation of the Somatic MAF file.  In a later release we will supply corrected BAM, VCF, and MAF files for these cases.  In the mean time, we advise you not to use any of the WXS files associated with these cases in the GDC Data Portal.  A list of these files can be found here. [Download list of affected files](DLBC_Affected_Files.txt). <!-- Data-871-->
    * TCGA-FF-8062
    * TCGA-FM-8000
    * TCGA-G8-6324
    * TCGA-G8-6325
    * TCGA-G8-6326
    * TCGA-G8-6906
    * TCGA-G8-6907
    * TCGA-G8-6909
    * TCGA-G8-6914
    * TCGA-GR-7351
    * TCGA-GR-7353
* Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* No data from TARGET-MDLS is available.
* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
* SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
* Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
* SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
* Tumor grade property is not populated <!--SV-585-->
* Progression_or_recurrence property is not populated <!--SV-584-->



## Data Release 10.0 <!--REQ-318-->

* __GDC Product__: Data
* __Release Date__: December 21, 2017

### New updates

1. New TARGET files for all projects <!--DAT-1135-->
2. TARGET updates for clinical and biospecimen data <!--DAT-1125, DAT-1153-->
3. Replace corrupted .bai files <!--TT-109-->
4. Update TCGA and TARGET MAF files to include VarScan2 indels and more information in all_effects column <!--DAT-987-->
5. Update VarScan VCF files <!--DAT-1020-->

Updated files for this release are listed [here](DR10.0_files_swap.txt.gz).
A complete list of files for DR10.0 are listed for the GDC Data Portal [here](gdc_manifest_20171221_data_release_10.0_active.txt.gz) and the GDC Legacy Archive [here](gdc_manifest_20171221_data_release_10.0_legacy.txt.gz).

### Bugs Fixed Since Last Release

None

### Known Issues and Workarounds

* The raw and annotated VarScan VCF files for aliquot `TCGA-VR-A8ET-01A-11D-A403-09` were not replaced in DR10.0 and thus do not contain indels.  However, the indels from this aliquot can be found in the MAF files and are displayed in the Exploration section in the Data Portal.  These VCFs files will be replaced in a later release.
<!--SV-950-->
* There are 5051 TARGET files for which `experimental_strategy`, `data_format`, `platform`, and `data_subtype` are blank <!--SV-944-->
* There are two cases with identical submitter_id `TARGET-10-PARUYU` <!--SV-940-->
* TARGET-MDLS cases do not have disease_type or primary_site populated <!--SV-939-->
* Some TARGET cases are missing `days_to_last_follow_up` <!--SV-934-->
* Some TARGET cases are missing `age_at_diagnosis` <!--SV-933-->
* Some TARGET files are not connected to all related aliquots <!--SV-929-->
* miRNA alignments include QC failed reads.  
* Samples of TARGET sample_type `Recurrent Blood Derived Cancer - Bone Marrow` are mislabeled as `Recurrent Blood Derived Cancer - Peripheral Blood`.  A workaround is to look at the sample barcode, which is -04 for `Recurrent Blood Derived Cancer - Bone Marrow`. (e.g. `TARGET-20-PAMYAS-04A-03R`) <!--SV-918-->
* FM-AD clinical and biospecimen supplement files have incorrect data format.  They are listed as XLSX, but are in fact TSV files. <!--DAT-1123-->
* Mutation frequency may be underestimated when using MAF files for genes that overlap other genes.  This is because MAF files only record one gene per variant.
* Most intronic mutations are removed for MAF generation.  However, validated variants may rescue these in some cases.  Therefore intronic mutations in MAF files are not representative of those called by mutation callers.
* The latest TARGET data is not yet available at the GDC.  For the complete and latest data, please see the [TARGET Data Matrix](https://ocg.cancer.gov/programs/target/data-matrix).  Data that is not present or is not the most up to date includes:
    *  All microarray data and metadata
    *  All sequencing analyzed data and metadata
    *  1180 of 12063 sequencing runs of raw data
* Demographic information for some TARGET patients is incorrect.  The correct information can be found in the associated clinical supplement file.  Impacted patients are TARGET-50-PAJNUS. <!--SV-710-->
* There are 11 cases in project TCGA-DLBC that are known to have incorrect WXS data in the GDC Data Portal.  Impacted cases are listed below.  This affects the BAMs and VCFs associated with these cases in the GDC Data Portal.  Corrected BAMs can be found in the GDC Legacy Archive.  Variants from affected aliquots appear in the protected MAFs with GDC_FILTER=ContEst to indicate a sample contamination problem, but are removed during the generation of the Somatic MAF file.  In a later release we will supply corrected BAM, VCF, and MAF files for these cases.  In the mean time, we advise you not to use any of the WXS files associated with these cases in the GDC Data Portal.  A list of these files can be found here. [Download list of affected files](DLBC_Affected_Files.txt). <!-- Data-871-->
    * TCGA-FF-8062
    * TCGA-FM-8000
    * TCGA-G8-6324
    * TCGA-G8-6325
    * TCGA-G8-6326
    * TCGA-G8-6906
    * TCGA-G8-6907
    * TCGA-G8-6909
    * TCGA-G8-6914
    * TCGA-GR-7351
    * TCGA-GR-7353
* Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* No data from TARGET-MDLS is available.
* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
* SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
* Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
* SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
* Tumor grade property is not populated <!--SV-585-->
* Progression_or_recurrence property is not populated <!--SV-584-->




## Data Release 9.0 <!--REQ-317-->

* __GDC Product__: Data
* __Release Date__: October 24, 2017

### New updates

1. Foundation Medicine Data Release
* This includes controlled-access VCF and MAF files as well as clinical and biospecimen supplements and metadata.
* Original Foundation Medicine supplied data can be found on the [Foundation Medicine Project Page](https://gdc.cancer.gov/about-gdc/contributed-genomic-data-cancer-research/foundation-medicine/foundation-medicine).
2.  Updated RNA-Seq data for TARGET NBL  
* Includes new BAM and count files

Updated files for this release are listed [here](DR9.0_files_swap.txt.gz).
A complete list of files for DR9.0 are listed [here](gdc_manifest_20171024_data_release_9.0_active.txt.gz).

### Bugs Fixed Since Last Release

None

### Known Issues and Workarounds

* miRNA alignments include QC failed reads.  
* Samples of TARGET sample_type `Recurrent Blood Derived Cancer - Bone Marrow` are mislabeled as `Recurrent Blood Derived Cancer - Peripheral Blood`.  A workaround is to look at the sample barcode, which is -04 for `Recurrent Blood Derived Cancer - Bone Marrow`. (e.g. `TARGET-20-PAMYAS-04A-03R`) <!--SV-918-->
* FM-AD clinical and biospecimen supplement files have incorrect data format.  They are listed as XLSX, but are in fact TSV files. <!--DAT-1123-->
* Mutation frequency may be underestimated when using MAF files for genes that overlap other genes.  This is because MAF files only record one gene per variant.
* Most intronic mutations are removed for MAF generation.  However, validated variants may rescue these in some cases.  Therefore intronic mutations in MAF files are not representative of those called by mutation callers.
* The latest TARGET data is not yet available at the GDC.  For the complete and latest data, please see the [TARGET Data Matrix](https://ocg.cancer.gov/programs/target/data-matrix).  Data that is not present or is not the most up to date includes:
    *  All microarray data and metadata
    *  All sequencing analyzed data and metadata
    *  1180 of 12063 sequencing runs of raw data
* Demographic information for some TARGET patients is incorrect.  The correct information can be found in the associated clinical supplement file.  Impacted patients are TARGET-50-PAJNUS. <!--SV-710-->
* There are 11 cases in project TCGA-DLBC that are known to have incorrect WXS data in the GDC Data Portal.  Impacted cases are listed below.  This affects the BAMs and VCFs associated with these cases in the GDC Data Portal.  Corrected BAMs can be found in the GDC Legacy Archive.  Variants from affected aliquots appear in the protected MAFs with GDC_FILTER=ContEst to indicate a sample contamination problem, but are removed during the generation of the Somatic MAF file.  In a later release we will supply corrected BAM, VCF, and MAF files for these cases.  In the mean time, we advise you not to use any of the WXS files associated with these cases in the GDC Data Portal.  A list of these files can be found here. [Download list of affected files](DLBC_Affected_Files.txt). <!-- Data-871-->
    * TCGA-FF-8062
    * TCGA-FM-8000
    * TCGA-G8-6324
    * TCGA-G8-6325
    * TCGA-G8-6326
    * TCGA-G8-6906
    * TCGA-G8-6907
    * TCGA-G8-6909
    * TCGA-G8-6914
    * TCGA-GR-7351
    * TCGA-GR-7353
* Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* No data from TARGET-MDLS is available.
* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
* SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
* Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
* SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
* Tumor grade property is not populated <!--SV-585-->
* Progression_or_recurrence property is not populated <!--SV-584-->

## Data Release 8.0

* __GDC Product__: Data
* __Release Date__: August 22, 2017

### New updates

1. Released updated miRNA quantification files to address double counting of some normalized counts described in DR7.0 release notes. <!--DAT-1, DAT-988-->

Updated files for this release are listed [here](DR8.0_files_swap.txt).
A Complete list of files for DR8.0 are listed [here](gdc_manifest_20170822_data_release_8.0_active.txt.gz).

### Bugs Fixed Since Last Release

None

### Known Issues and Workarounds

* TARGET-NBL RNA-Seq files were run as single ended even though they are derived from paired-end data.  These files will be rerun through the GDC RNA-Seq pipelines in a later release.  Impacted files can be found [here](https://portal.gdc.cancer.gov/repository?facetTab=files&files_offset=20&filters=~%28op~%27and~content~%28~%28op~%27in~content~%28field~%27cases.project.program.name~value~%28~%27TARGET%29%29%29~%28op~%27in~content~%28field~%27cases.project.project_id~value~%28~%27TARGET-NBL%29%29%29~%28op~%27in~content~%28field~%27files.data_format~value~%28~%27BAM%29%29%29~%28op~%27in~content~%28field~%27files.experimental_strategy~value~%28~%27RNA-Seq%29%29%29%29%29&searchTableTab=files).  Downstream count files are also affected.  Users may access original FASTQ files in the GDC Legacy Archive, which are not impacted by this issue.
* Mutation frequency may be underestimated when using MAF files for genes that overlap other genes.  This is because MAF files only record one gene per variant.
* Most intronic mutations are removed for MAF generation.  However, validated variants may rescue these in some cases.  Therefore intronic mutations in MAF files are not representative of those called by mutation callers.
* The latest TARGET data is not yet available at the GDC.  For the complete and latest data, please see the [TARGET Data Matrix](https://ocg.cancer.gov/programs/target/data-matrix).  Data that is not present or is not the most up to date includes:
    *  All microarray data and metadata
    *  All sequencing analyzed data and metadata
    *  1180 of 12063 sequencing runs of raw data
* Demographic information for some TARGET patients is incorrect.  The correct information can be found in the associated clinical supplement file.  Impacted patients are TARGET-50-PAJNUS. <!--SV-710-->
* There are 11 cases in project TCGA-DLBC that are known to have incorrect WXS data in the GDC Data Portal.  Impacted cases are listed below.  This affects the BAMs and VCFs associated with these cases in the GDC Data Portal.  Corrected BAMs can be found in the GDC Legacy Archive.  Variants from affected aliquots appear in the protected MAFs with GDC_FILTER=ContEst to indicate a sample contamination problem, but are removed during the generation of the Somatic MAF file.  In a later release we will supply corrected BAM, VCF, and MAF files for these cases.  In the mean time, we advise you not to use any of the WXS files associated with these cases in the GDC Data Portal.  A list of these files can be found here. [Download list of affected files](DLBC_Affected_Files.txt). <!-- Data-871-->
    * TCGA-FF-8062
    * TCGA-FM-8000
    * TCGA-G8-6324
    * TCGA-G8-6325
    * TCGA-G8-6326
    * TCGA-G8-6906
    * TCGA-G8-6907
    * TCGA-G8-6909
    * TCGA-G8-6914
    * TCGA-GR-7351
    * TCGA-GR-7353
* Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* No data from TARGET-MDLS is available.
* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
* SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
* Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
* SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
* Tumor grade property is not populated <!--SV-585-->
* Progression_or_recurrence property is not populated <!--SV-584-->




## Data Release 7.0

* __GDC Product__: Data
* __Release Date__: June 29, 2017

### New updates

1.  Updated public Mutation Annotation Format (MAF) files are now available. Updates include filtering to remove variants impacted by OxoG artifacts and those impacted by strand bias. <!--AP-29-->
2.  Protected MAF files are updated to include flags for OxoG and strand bias. <!--AP-29-->
2.  Annotated VCFs are updated to include flags for OxoG artifacts and strand bias. <!--DAT-874-->

Updated files for this release are listed [here](DR7.0_files_swap.txt).
A Complete list of files for DR7.0 are listed [here](gdc_manifest_20170629_data_release_7.0.txt.gz)

### Bugs Fixed Since Last Release

None

### Known Issues and Workarounds

* TARGET-NBL RNA-Seq files were run as single ended even though they are derived from paired-end data.  These files will be rerun through the GDC RNA-Seq pipelines in a later release.  Impacted files can be found [here](https://portal.gdc.cancer.gov/repository?facetTab=files&files_offset=20&filters=~%28op~%27and~content~%28~%28op~%27in~content~%28field~%27cases.project.program.name~value~%28~%27TARGET%29%29%29~%28op~%27in~content~%28field~%27cases.project.project_id~value~%28~%27TARGET-NBL%29%29%29~%28op~%27in~content~%28field~%27files.data_format~value~%28~%27BAM%29%29%29~%28op~%27in~content~%28field~%27files.experimental_strategy~value~%28~%27RNA-Seq%29%29%29%29%29&searchTableTab=files).  Downstream count files are also affected.  Users may access original FASTQ files in the GDC Legacy Archive, which are not impacted by this issue.
* Reads that are mapped to multiple genomic locations are double counted in some of the GDC miRNA results.  The GDC will release updated files correcting the issue in an upcoming release.<!--SV-745-->  The specific impacts are described further below:
    *  Isoform Expression Quantification files
        *  Raw reads counts are accurate
        *  Normalized counts are proportionally skewed (r^2=1.0)
    *  miRNA Expression Quantification files
        *  A small proportion of miRNA counts are overestimated (mean r^2=0.9999)
        *  Normalized counts are proportionally skewed (mean r^2=0.9999)
    * miRNA BAM files
        * no impact
* Mutation frequency may be underestimated when using MAF files for genes that overlap other genes.  This is because MAF files only record one gene per variant.
* Most intronic mutations are removed for MAF generation.  However, validated variants may rescue these in some cases.  Therefore intronic mutations in MAF files are not representative of those called by mutation callers.
* The latest TARGET data is not yet available at the GDC.  For the complete and latest data, please see the [TARGET Data Matrix](https://ocg.cancer.gov/programs/target/data-matrix).  Data that is not present or is not the most up to date includes:
    *  All microarray data and metadata
    *  All sequencing analyzed data and metadata
    *  1180 of 12063 sequencing runs of raw data
* Demographic information for some TARGET patients is incorrect.  The correct information can be found in the associated clinical supplement file.  Impacted patients are TARGET-50-PAJNUS. <!--SV-710-->
* There are 11 cases in project TCGA-DLBC that are known to have incorrect WXS data in the GDC Data Portal.  Impacted cases are listed below.  This affects the BAMs and VCFs associated with these cases in the GDC Data Portal.  Corrected BAMs can be found in the GDC Legacy Archive.  Variants from affected aliquots appear in the protected MAFs with GDC_FILTER=ContEst to indicate a sample contamination problem, but are removed during the generation of the Somatic MAF file.  In a later release we will supply corrected BAM, VCF, and MAF files for these cases.  In the mean time, we advise you not to use any of the WXS files associated with these cases in the GDC Data Portal.  A list of these files can be found here. [Download list of affected files](DLBC_Affected_Files.txt). <!-- Data-871-->
    * TCGA-FF-8062
    * TCGA-FM-8000
    * TCGA-G8-6324
    * TCGA-G8-6325
    * TCGA-G8-6326
    * TCGA-G8-6906
    * TCGA-G8-6907
    * TCGA-G8-6909
    * TCGA-G8-6914
    * TCGA-GR-7351
    * TCGA-GR-7353
* Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* No data from TARGET-MLDS is available.
* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
* SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
* Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
* SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
* Tumor grade property is not populated <!--SV-585-->
* Progression_or_recurrence property is not populated <!--SV-584-->




## Data Release 6.0

* __GDC Product__: Data
* __Release Date__: May 9, 2017

### New updates

1.  GDC updated public Mutation Annotation Format (MAF) files are now available. Updates include leveraging the MC3 variant filtering strategy, which results in more variants being recovered relative to the previous version. A detailed description of the new format can be found [here](../File_Formats/MAF_Format/). <!--DAT-572-->
2.  Protected MAFs are updated to include additional variant annotation information <!--DAT-572-->
3.  Some MuTect2 VCFs updated to include dbSNP and COSMIC annotations found in other VCFs <!--TT-21-->

Updated files for this release are listed [here](DR6.0_files_swap.txt).

### Bugs Fixed Since Last Release

None

### Known Issues and Workarounds

* There are 11 cases in project TCGA-DLBC that are known to have incorrect WXS data in the GDC Data Portal.  Impacted cases are listed below.  This affects the BAMs and VCFs associated with these cases in the GDC Data Portal.  Corrected BAMs can be found in the GDC Legacy Archive.  Variants from affected aliquots appear in the protected MAFs with GDC_FILTER=ContEst to indicate a sample contamination problem, but are removed during the generation of the Somatic MAF file.  In a later release we will supply corrected BAM, VCF, and MAF files for these cases.  In the mean time, we advise you not to use any of the WXS files associated with these cases in the GDC Data Portal.  A list of these files can be found here. [Download list of affected files](DLBC_Affected_Files.txt). <!-- Data-871-->
    * TCGA-FF-8062
    * TCGA-FM-8000
    * TCGA-G8-6324
    * TCGA-G8-6325
    * TCGA-G8-6326
    * TCGA-G8-6906
    * TCGA-G8-6907
    * TCGA-G8-6909
    * TCGA-G8-6914
    * TCGA-GR-7351
    * TCGA-GR-7353
* Variants found in VCF and MAF files may contain OxoG artifacts, which are produced during library preparation and may result in the apparent substitutions of C to A or G to T in certain sequence contexts.  In the future we will plan to label potential oxoG artifacts in the MAF files.
* Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
* Some validated somatic mutations may not be present in open-access MAF files.  Please review the protected MAF files in the GDC Data Portal if you are unable to find your mutation in the open-access files.
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* No data from TARGET-MLDS is available.
* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
* SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
* Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
* SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
* Tumor grade property is not populated <!--SV-585-->
* Progression_or_recurrence property is not populated <!--SV-584-->


Details are provided in [Data Release Manifest](Manifests/GDC_Data_v3_release_notes_manifest.txt)
<br>

## Data Release 5.0

* __GDC Product__: Data
* __Release Date__: March 16, 2017

### New updates

1.  Additional annotations from TCGA DCC are available <!--DAT-52-->
    * Complete list of updated TCGA files is found [here](DR5.0_CHANGES_TCGA.xlsx)
2.  Clinical data added for TARGET ALL P1 and P2 <!--DAT-197-->
3.  Pathology reports now have submitter IDs as assigned by the BCR <!--DAT-81-->
4.  TARGET Data refresh
    * Most recent biospecimen and clinical information from the TARGET DCC. New imported files are listed [here](DR5.0_changes_TARGET.xlsx)
    * Updated indexed biospecimen and clinical metadata
    * Updated SRA XMLs files
    * Does not include updates to TARGET NBL <!--SV-585-->

### Bugs Fixed Since Last Release

1. Missing cases from TCGA-LAML were added to Legacy Archive <!--DAT-272,DAT-324-->
2. Biotab files are now linked to Projects and Cases in Legacy Archive <!--SV-303,DAT-28-->

### Known Issues and Workarounds

* Some TCGA annotations are unavailable in the Legacy Archive or Data Portal<!--DAT-52-->. These annotations can be found [here](tcga-annotations-unavailable-20170315.json).
* Some validated somatic mutations may not be present in open-access MAF files.  When creating open-access MAF files from the protected versions we are extremely conservative in removing potential germline variants.  Our approach is to remove all mutations that are present in dbSNP.  In a subsequent release we will provide updated open-access MAF files, which preserve variants found in MC3 or a TCGA validation study.  Please review the protected MAF files in the GDC Data Portal if you are unable to find your mutation in the open-access files.
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* MAF Column #109 "FILTER" entries are separated by both commas and semi-colons. <!-- PGDC-2589 -->
* TARGET-AML is undergoing reorganization.  Pending reorganization, cases from this projects may not contain many clinical, biospecimen, or genomic data files.
* No data from TARGET-MLDS is available.
* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
* SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
* Two biotab files are not linked to Project or Case in the Legacy Archive <!--SV-535, DAT-493-->
* SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->
* Tumor grade property is not populated <!--SV-585-->
* Progression_or_recurrence property is not populated <!--SV-584-->


Details are provided in [Data Release Manifest](Manifests/GDC_Data_v3_release_notes_manifest.txt)
<br>


## Data Release 4.0

* __GDC Product__: Data
* __Release Date__: October 31, 2016

### New updates

1.  TARGET ALL P1 and P2 biospecimen and molecular data are now available in the Legacy Archive.  Clinical data will be available in a later release. <!-- Dat-185, Dat-194-->
2.  Methylation data from 27k/450k Arrays has been lifted over to hg38 and is now available in the GDC Data Portal <!-- Dat-109 -->
3.  Public MAF files are now available for VarScan2, MuSE, and SomaticSniper.  MuTect2 MAFs were made available in a previous release. <!--DAT-235-->
4.  Updated VCFs and MAF files are available for MuTect2 pipeline to compensate for WGA-related false positive indels.  See additional information on that change [here](https://gdc.cancer.gov/content/mutect2-insertion-artifacts). A listing of replaced files is provided [here](Manifests/GDC_Data_v4_mapping_of_replaced_Mutect2_MAF_and_VCF_files.zip). <!-- Dat-145, Dat-260 -->
5.  Added submitter_id for Pathology Reports in Legacy Archive <!--DAT-81-->

### Bugs Fixed Since Last Release

* None

### Known Issues and Workarounds

* Some validated somatic mutations may not be present in open-access MAF files.  When creating open-access MAF files from the protected versions we are extremely conservative in removing potential germline variants.  Our approach is to remove all mutations that are present in dbSNP.  In a subsequent release we will provide updated open-access MAF files, which preserve variants found in COSMIC or a TCGA validation study.  Please review the protected MAF files in the GDC Data Portal if you are unable to find your mutation in the open-access files.
* Public MAF files for different variant calling pipelines but the same project may contain different numbers of samples.  Samples are omitted from the public MAF files if they have no PASS variants, which can lead to this apparent discrepancy.
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* MAF Column #109 "FILTER" entries are separated by both commas and semi-colons. <!-- PGDC-2589 -->
* TARGET-AML is undergoing reorganization.  Pending reorganization, cases from this projects may not contain many clinical, biospecimen, or genomic data files.
* No data from TARGET-MLDS is available.
* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
* SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
* There are 200 cases from TCGA-LAML that do not appear in the Legacy Archive <!--SV-327-->
* Biotab files are not linked to Project or Case in the Legacy Archive <!--SV-303-->
* SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->


Details are provided in [Data Release Manifest](Manifests/GDC_Data_v5_release_notes_manifest.txt)
<br>



## Data Release 3.0

* __GDC Product__: Data
* __Release Date__: September 16, 2016

### New updates

1.  CCLE data now available (in the Legacy Archive only)
2.  BMI calculation is corrected
3.  Slide is now categorized as a Biospecimen entity

### Bugs Fixed Since Last Release

* BMI calculation is corrected

### Known Issues and Workarounds

* Insertions called for tumor samples that underwent whole genome amplification may be of lower quality.  Whether a sample underwent this process can be found in the analyte_type property within analyte and aliquot. TCGA analyte type can be also identified in the 20th character of TCGA barcode, at which "W" corresponds to WGA.<!--BINF-6-->
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* Public MAFs (those with germline variants removed) are only available for MuTect2 pipeline.  MAFs for other pipelines are forthcoming.  
* MAF Column #109 "FILTER" entries are separated by both commas and semi-colons. <!-- PGDC-2589 -->
* TARGET-AML and TARGET-ALL projects are undergoing reorganization.  Pending reorganization, cases from these projects may not contain many clinical, biospecimen, or genomic data files.
* No data from TARGET-PPTP is available.
* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
* SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
* There are 200 cases from TCGA-LAML that do not appear in the Legacy Archive <!--SV-327-->
* Biotab files are not linked to Project or Case in the Legacy Archive <!--SV-303-->
* SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->


Details are provided in [Data Release Manifest](Manifests/GDC_Data_v3_release_notes_manifest.txt)
<br>


## Data Release 2.0

* __GDC Product__: Data
* __Release Date__: August 9, 2016


### New updates

1.  Additional data, previously available via CGHub and the TCGA DCC, is now available in the GDC
2.  Better linking between files and their associated projects and cases in the Legacy Archive
3.  MAF files are now available in the GDC Data Portal

### Known Issues and Workarounds

* Insertions called for tumor samples that underwent whole genome amplification may be of lower quality.  These are present in VCF and MAF files produced by the MuTect2 variant calling pipeline.  This information can be found in the analyte_type property within analyte and aliquot. TCGA analyte type can be also identified in the 20th character of TCGA barcode, at which "W" corresponds to WGA.<!--BINF-6-->
* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* Public MAFs (those with germline variants removed) are only available for MuTect2 pipeline.  MAFs for other pipelines are forthcoming.  
* MAF Column #109 "FILTER" entries are separated by both commas and semi-colons. <!-- PGDC-2589 -->
* TARGET-AML and TARGET-ALL projects are undergoing reorganization.  Pending reorganization, cases from these projects may not contain many clinical, biospecimen, or genomic data files.
* No data from TARGET-PPTP is available.
* Slide barcodes (`submitter_id` values for Slide entities in the Legacy Archive) are not available <!-- DAT-10 -->
* SDF Files are not linked to Project or Case in the Legacy Archive <!--SV-332-->
* There are 200 cases from TCGA-LAML that do not appear in the Legacy Archive <!--SV-327-->
* Biotab files are not linked to Project or Case in the Legacy Archive <!--SV-303-->
* SDRF files are not linked to Project or Case in the Legacy Archive <!--SV-288-->
* Portion "weight" property is incorrectly described in the Data Dictionary as the weight of the patient in kg, should be described as the weight of the portion in mg <!--SV-391-->

Details are provided in [Data Release Manifest](Manifests/GDC_Data_v2_release_notes_manifest.txt)
<br>

## Initial Data Release (1.0)

* __GDC Product__: Data
* __Release Date__: June 6, 2016

### Available Program Data

* The Cancer Genome Atlas (TCGA)
* Therapeutically Applicable Research To Generate Effective Treatments (TARGET)

### Available Harmonized Data

* WXS
     * Co-cleaned BAM files aligned to GRCh38 using BWA
* mRNA-Seq
    * BAM files aligned to GRCh38 using STAR 2-pass strategy
    * Expression quantification using HTSeq
* miRNA-Seq
    * BAM files aligned to GRCh38 using BWA aln
    * Expression quantification using BCCA miRNA Profiling Pipeline*
* Genotyping Array
    * CNV segmentation data


### Known Issues and Workarounds

* BAM files produced by the GDC RNA-Seq Alignment workflow will currently fail validation using the Picard ValidateSamFiles tool.  This is caused by STAR2 not recording mate mapping information for unmapped reads, which are retained in our BAM files.  Importantly, all affected BAM files are known to behave normally in downstream workflows including expression quantification.
* All legacy files for TCGA are available in the GDC Legacy Archive, but not always linked back to cases depending on available metadata.
* Public MAFs (those with germline variants removed) are only available for MuTect2 pipeline.  MAFs for other pipelines are forthcoming.  
* TARGET-AML and TARGET-ALL projects are undergoing reorganization.  Pending reorganization, cases from these projects may not contain many clinical, biospecimen, or genomic data files.
* No data from TARGET-PPTP is available.
*	Legacy data not available in harmonized form:
    *	Annotated VCF files from TARGET, anticipated in future data release
    * TCGA data that failed harmonization or QC or have been newly updated in CGHub: ~1.0% of WXS aliquots, ~1.6% of RNA-Seq aliquots
    * TARGET data that failed harmonization or QC, have been newly updated in CGHub, or whose project names are undergoing reorganization: ~76% of WXS aliquots, ~49% of RNA-Seq aliquots, ~57% of miRNA-Seq.
* MAF Column #109 "FILTER" entries are separated by both commas and semi-colons. <!-- PGDC-2589 -->
*	MAFs are not yet available for query or search in the GDC Data Portal or API.  You may download these files using the following manifests, which can be passed directly to the Data Transfer Tool.  Links for the open-access TCGA MAFs are provided below for downloading individual files.
    * [Open-access MAFs manifest](Manifests/GDC_open_MAFs_manifest.txt)
    * [Controlled-access MAFs manifest](Manifests/GDC_controlled_MAFs_manifest.txt)

Details are provided in [Data Release Manifest](Manifests/GDC_Data_v1_release_notes_manifest.txt)

### Download Open-access MAF files

* Please note that these links no longer point to files and will be updated in the future.

<a href="https://api.gdc.cancer.gov/data/abbe72a5-cb39-48e4-8df5-5fd2349f2bb2">TCGA.ACC.mutect.somatic.maf.gz</a><br>
<a href="https://api.gdc.cancer.gov/data/cb97ef6e-7a13-4d42-81b4-1510c00d6373">TCGA.BLCA.mutect.somatic.maf.gz</a><br>
<a href="https://api.gdc.cancer.gov/data/91ae5ca9-55c2-4c9c-929e-8638444dc7b5">TCGA.BRCA.mutect.somatic.maf.gz</a><br>
<a href="https://api.gdc.cancer.gov/data/fedc0e90-5070-41e4-993d-603b92f4ecfd">TCGA.CESC.mutect.somatic.maf.gz</a><br>
<a href="https://api.gdc.cancer.gov/data/f33fd38b-287c-4978-a1bb-a95bbfd4351a">TCGA.CHOL.mutect.somatic.maf.gz</a><br>
<a href="https://api.gdc.cancer.gov/data/bf1c53dc-79bb-43ae-88e4-23758853e5c6">TCGA.COAD.mutect.somatic.maf.gz</a><br>
<a href="https://api.gdc.cancer.gov/data/4835f959-6ab4-4ee8-901f-c92aaad4592d">TCGA.DLBC.mutect.somatic.maf.gz</a><br>
<a href="https://api.gdc.cancer.gov/data/150c6a9f-71cd-4710-9617-cd150498202e">TCGA.ESCA.mutect.somatic.maf.gz</a><br>
<a href="https://api.gdc.cancer.gov/data/4b176a7b-a5c3-457e-af95-992018b6f3d7">TCGA.GBM.mutect.somatic.maf.gz</a><br>
<a href="https://api.gdc.cancer.gov/data/64683606-b957-4478-a7d5-673de68b0341">TCGA.HNSC.mutect.somatic.maf.gz</a><br>
<a href="https://api.gdc.cancer.gov/data/a8dc2dd2-74b3-4035-9551-c0ae3f76293e">TCGA.KICH.mutect.somatic.maf.gz</a><br>
<a href="https://api.gdc.cancer.gov/data/0cf0f121-4f10-436e-bfc3-8fcfd5f78d0d">TCGA.KIRC.mutect.somatic.maf.gz</a><br>
<a href="https://api.gdc.cancer.gov/data/b1c13b93-7ad6-4d09-b613-84a7c55e61d9">TCGA.KIRP.mutect.somatic.maf.gz</a><br>
<a href="https://api.gdc.cancer.gov/data/bac5a617-5b1d-4c33-b9ac-b48bd7e4947a">TCGA.LAML.mutect.somatic.maf.gz</a><br>
<a href="https://api.gdc.cancer.gov/data/42ff7a98-5a9a-48ad-ad9d-d3a23c245296">TCGA.LGG.mutect.somatic.maf.gz</a><br>
<a href="https://api.gdc.cancer.gov/data/dc3f239a-ffc6-4e60-b5f5-9f365daaf60a">TCGA.LIHC.mutect.somatic.maf.gz</a><br>
<a href="https://api.gdc.cancer.gov/data/76b0eec4-bbb6-4340-972d-05a5aace63a4">TCGA.LUAD.mutect.somatic.maf.gz</a><br>
<a href="https://api.gdc.cancer.gov/data/846c4788-cda0-4f11-b240-a7ad977e032f">TCGA.LUSC.mutect.somatic.maf.gz</a><br>
<a href="https://api.gdc.cancer.gov/data/34a32349-bf87-4e96-86a5-ca23f0db475e">TCGA.MESO.mutect.somatic.maf.gz</a><br>
<a href="https://api.gdc.cancer.gov/data/bc5ee1aa-969d-472d-920b-0e654cc585fa">TCGA.OV.mutect.somatic.maf.gz</a><br>
<a href="https://api.gdc.cancer.gov/data/6f72cd49-6c0e-4409-8e94-26ea5d421bc8">TCGA.PAAD.mutect.somatic.maf.gz</a><br>
<a href="https://api.gdc.cancer.gov/data/e087cf64-b514-4e92-af9d-2b18341098d5">TCGA.PCPG.mutect.somatic.maf.gz</a><br>
<a href="https://api.gdc.cancer.gov/data/0d708001-437e-4b78-9ffa-bfafdfc10a28">TCGA.PRAD.mutect.somatic.maf.gz</a><br>
<a href="https://api.gdc.cancer.gov/data/221f1dec-b539-4345-b687-435659fc21af">TCGA.READ.mutect.somatic.maf.gz</a><br>
<a href="https://api.gdc.cancer.gov/data/84fad8c0-ac06-4181-92d9-0562392325ba">TCGA.SARC.mutect.somatic.maf.gz</a><br>
<a href="https://api.gdc.cancer.gov/data/4aee3e32-2802-4e1e-8577-d74b414f30f7">TCGA.SKCM.mutect.somatic.maf.gz</a><br>
<a href="https://api.gdc.cancer.gov/data/98d14107-5fb5-49bd-ac38-a52178838d6c">TCGA.STAD.mutect.somatic.maf.gz</a><br>
<a href="https://api.gdc.cancer.gov/data/c3155e81-df07-49c6-8502-ef6ebc60812e">TCGA.TGCT.mutect.somatic.maf.gz</a><br>
<a href="https://api.gdc.cancer.gov/data/3207f158-6f79-4636-81a1-ce1b30157925">TCGA.THCA.mutect.somatic.maf.gz</a><br>
<a href="https://api.gdc.cancer.gov/data/ea93d4df-5e76-484a-b7b8-93900ea1d61c">TCGA.THYM.mutect.somatic.maf.gz</a><br>
<a href="https://api.gdc.cancer.gov/data/a142d8b8-7741-4869-9ca4-0025890eee18">TCGA.UCEC.mutect.somatic.maf.gz</a><br>
<a href="https://api.gdc.cancer.gov/data/cfc3f9e4-ac0d-4133-b616-a7d2caf28e54">TCGA.UCS.mutect.somatic.maf.gz</a><br>
<a href="https://api.gdc.cancer.gov/data/a96185ef-5b9b-4f0e-b437-f6a0f4f0892b">TCGA.UVM.mutect.somatic.maf.gz</a><br>
