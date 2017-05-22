# Data Release Notes

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

* There are 11 cases in project TCGA-DLBC that are known to have incorrect WXS data in the GDC Data Portal.  Impacted cases are listed below.  This affects the BAMs and VCFs associated with these cases in the GDC Data Portal.  Corrected BAMs can be found in the GDC Legacy Archive.  Variants from affected aliquots appear in the protected MAFs with GDC_FILTER=ContEst to indicate sample contamination problem, but are removed during the generation of the Somatic MAF file.  In a later release we will supply corrected BAM, VCF, and MAF files for these cases.  In the mean time, we advise you not to use any of the WXS files associated with these cases in the GDC Data Portal ([DLBC Affected Files](DLBC_Affected_Files.txt)). <!-- Data-871-->
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
4.  Updated VCFs and MAF files are available for MuTect2 pipeline to compensate for WGA-related false positive indels.  See additional information on that change [here](https://gdc.cancer.gov/about-gdc/scientific-reports/known-mutect2-variant-artifacts). A listing of replaced files is provided [here](GDC_Data_v4_mapping_of_replaced_Mutect2_MAF_and_VCF_files.zip). <!-- Dat-145, Dat-260 -->
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
