# Data Release Notes

| Version | Date |
|---|---|
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
4.  Updated VCFs and MAF files are available for MuTect2 pipeline to compensate for WGA-related false positive indels.  See additional information on that change [here](https://gdc.cancer.gov/content/mutect2-insertion-artifacts). A listing of replaced files is provided [here](GDC_Data_v4_mapping_of_replaced_Mutect2_MAF_and_VCF_files.zip). <!-- Dat-145, Dat-260 -->
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
