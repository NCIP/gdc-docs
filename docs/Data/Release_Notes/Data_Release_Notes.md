# Data Release Notes






## Initial Data Release

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
* TARGET-AML and TARGET-ALL projects are undergoing reorganization.  Pending reorganization, cases from these projects may not contain many Clinical, Biospecimen, or genomic data files.
*	Legacy data not available in harmonized form:
    *	Annotated VCF files from TARGET, anticipated in future data release
    * TCGA data that failed harmonization or QC or have been newly updated in CGHub: ~1.0% of WXS aliquots, ~1.6% of RNA-Seq aliquots
    * All TARGET-PPTP data
    * TARGET data that failed harmonization or QC, have been newly updated in CGHub, or whose project names are undergoing reorganization: ~76% of WXS aliquots, ~49% of RNA-Seq aliquots, ~57% of miRNA-Seq.
*	MAFs are not yet available for query or search in the GDC Data Portal or API.  You may download these files using the following manifests, which can be passed directly to the Data Transfer Tool.  Links for the open-access MAFs are provided below for downloading invidual files.
    * [Open-access MAFs manifest](Manifests/GDC_open_MAFs_manifest.txt)
    * [Controlled-access MAFs manifest](Manifests/GDC_controlled_MAFs_manifest.txt)

Details are provided in [Data Release Manifest](Manifests/GDC_Data_v1_release_notes_manifest.txt)

### Download Open-access MAF files
<a href="https://gdc-api.nci.nih.gov/data/abbe72a5-cb39-48e4-8df5-5fd2349f2bb2">TCGA.ACC.mutect.somatic.maf.gz</a><br>
<a href="https://gdc-api.nci.nih.gov/data/cb97ef6e-7a13-4d42-81b4-1510c00d6373">TCGA.BLCA.mutect.somatic.maf.gz</a><br>
<a href="https://gdc-api.nci.nih.gov/data/91ae5ca9-55c2-4c9c-929e-8638444dc7b5">TCGA.BRCA.mutect.somatic.maf.gz</a><br>
<a href="https://gdc-api.nci.nih.gov/data/fedc0e90-5070-41e4-993d-603b92f4ecfd">TCGA.CESC.mutect.somatic.maf.gz</a><br>
<a href="https://gdc-api.nci.nih.gov/data/f33fd38b-287c-4978-a1bb-a95bbfd4351a">TCGA.CHOL.mutect.somatic.maf.gz</a><br>
<a href="https://gdc-api.nci.nih.gov/data/bf1c53dc-79bb-43ae-88e4-23758853e5c6">TCGA.COAD.mutect.somatic.maf.gz</a><br>
<a href="https://gdc-api.nci.nih.gov/data/4835f959-6ab4-4ee8-901f-c92aaad4592d">TCGA.DLBC.mutect.somatic.maf.gz</a><br>
<a href="https://gdc-api.nci.nih.gov/data/150c6a9f-71cd-4710-9617-cd150498202e">TCGA.ESCA.mutect.somatic.maf.gz</a><br>
<a href="https://gdc-api.nci.nih.gov/data/4b176a7b-a5c3-457e-af95-992018b6f3d7">TCGA.GBM.mutect.somatic.maf.gz</a><br>
<a href="https://gdc-api.nci.nih.gov/data/64683606-b957-4478-a7d5-673de68b0341">TCGA.HNSC.mutect.somatic.maf.gz</a><br>
<a href="https://gdc-api.nci.nih.gov/data/a8dc2dd2-74b3-4035-9551-c0ae3f76293e">TCGA.KICH.mutect.somatic.maf.gz</a><br>
<a href="https://gdc-api.nci.nih.gov/data/0cf0f121-4f10-436e-bfc3-8fcfd5f78d0d">TCGA.KIRC.mutect.somatic.maf.gz</a><br>
<a href="https://gdc-api.nci.nih.gov/data/b1c13b93-7ad6-4d09-b613-84a7c55e61d9">TCGA.KIRP.mutect.somatic.maf.gz</a><br>
<a href="https://gdc-api.nci.nih.gov/data/bac5a617-5b1d-4c33-b9ac-b48bd7e4947a">TCGA.LAML.mutect.somatic.maf.gz</a><br>
<a href="https://gdc-api.nci.nih.gov/data/42ff7a98-5a9a-48ad-ad9d-d3a23c245296">TCGA.LGG.mutect.somatic.maf.gz</a><br>
<a href="https://gdc-api.nci.nih.gov/data/dc3f239a-ffc6-4e60-b5f5-9f365daaf60a">TCGA.LIHC.mutect.somatic.maf.gz</a><br>
<a href="https://gdc-api.nci.nih.gov/data/76b0eec4-bbb6-4340-972d-05a5aace63a4">TCGA.LUAD.mutect.somatic.maf.gz</a><br>
<a href="https://gdc-api.nci.nih.gov/data/846c4788-cda0-4f11-b240-a7ad977e032f">TCGA.LUSC.mutect.somatic.maf.gz</a><br>
<a href="https://gdc-api.nci.nih.gov/data/34a32349-bf87-4e96-86a5-ca23f0db475e">TCGA.MESO.mutect.somatic.maf.gz</a><br>
<a href="https://gdc-api.nci.nih.gov/data/bc5ee1aa-969d-472d-920b-0e654cc585fa">TCGA.OV.mutect.somatic.maf.gz</a><br>
<a href="https://gdc-api.nci.nih.gov/data/6f72cd49-6c0e-4409-8e94-26ea5d421bc8">TCGA.PAAD.mutect.somatic.maf.gz</a><br>
<a href="https://gdc-api.nci.nih.gov/data/e087cf64-b514-4e92-af9d-2b18341098d5">TCGA.PCPG.mutect.somatic.maf.gz</a><br>
<a href="https://gdc-api.nci.nih.gov/data/0d708001-437e-4b78-9ffa-bfafdfc10a28">TCGA.PRAD.mutect.somatic.maf.gz</a><br>
<a href="https://gdc-api.nci.nih.gov/data/221f1dec-b539-4345-b687-435659fc21af">TCGA.READ.mutect.somatic.maf.gz</a><br>
<a href="https://gdc-api.nci.nih.gov/data/84fad8c0-ac06-4181-92d9-0562392325ba">TCGA.SARC.mutect.somatic.maf.gz</a><br>
<a href="https://gdc-api.nci.nih.gov/data/4aee3e32-2802-4e1e-8577-d74b414f30f7">TCGA.SKCM.mutect.somatic.maf.gz</a><br>
<a href="https://gdc-api.nci.nih.gov/data/98d14107-5fb5-49bd-ac38-a52178838d6c">TCGA.STAD.mutect.somatic.maf.gz</a><br>
<a href="https://gdc-api.nci.nih.gov/data/c3155e81-df07-49c6-8502-ef6ebc60812e">TCGA.TGCT.mutect.somatic.maf.gz</a><br>
<a href="https://gdc-api.nci.nih.gov/data/3207f158-6f79-4636-81a1-ce1b30157925">TCGA.THCA.mutect.somatic.maf.gz</a><br>
<a href="https://gdc-api.nci.nih.gov/data/ea93d4df-5e76-484a-b7b8-93900ea1d61c">TCGA.THYM.mutect.somatic.maf.gz</a><br>
<a href="https://gdc-api.nci.nih.gov/data/a142d8b8-7741-4869-9ca4-0025890eee18">TCGA.UCEC.mutect.somatic.maf.gz</a><br>
<a href="https://gdc-api.nci.nih.gov/data/cfc3f9e4-ac0d-4133-b616-a7d2caf28e54">TCGA.UCS.mutect.somatic.maf.gz</a><br>
<a href="https://gdc-api.nci.nih.gov/data/a96185ef-5b9b-4f0e-b437-f6a0f4f0892b">TCGA.UVM.mutect.somatic.maf.gz</a><br>
