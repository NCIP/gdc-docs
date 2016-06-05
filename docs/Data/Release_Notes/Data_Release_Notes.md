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
*	Legacy data not available in harmonized form:
    *	Annotated VCF files from TARGET, anticipated in future data release
    * TCGA data that failed harmonization or QC or have been newly updated in CGHub: ~1.0% of WXS aliquots, ~1.6% of RNA-Seq aliquots
    * TARGET data that failed harmonization or QC, have been newly updated in CGHub, or whose project names are subject to change: ~76% of WXS aliquots, ~49% of RNA-Seq aliquots, ~57% of miRNA-Seq.
*	MAFs are not yet available for query or search in the GDC Data Portal or API.  You may download these files using the following manifests, which can be passed directly to the Data Transfer Tool.  Individual UUIDs can also be passed directly to the Data Transfer Tool or to the API (eg. https://gdc-api.nci.nih.gov/data/abbe72a5-cb39-48e4-8df5-5fd2349f2bb2).
    * [Open-access MAFs manifest](Manifests/GDC_open_MAFs_manifest.txt)
    * [Controlled-access MAFs manifest](Manifests/GDC_controlled_MAFs_manifest.txt)
* TARGET-AML and TARGET-ALL projects are undergoing reorganization.  Pending reorganization, cases from these projects may not contain many Clinical, Biospecimen, or genomic data files.

Details are provided in [Data Release Manifest](Manifests/GDC_Data_v1_release_notes_manifest.txt)
