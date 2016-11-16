# cBioPortal Release Notes

## Description

cBioPortal for Cancer Genomics provides visualization, analysis, and download of large-scale cancer genomics data sets.  The GDC hosts an instance of cBioPortal that displays harmonized data from TCGA powered by data from the GDC open-access MAF files. Visualizations related to RPPA, RNA-Seq, copy number variation, clinical data, as well as many external databases are not supported in the GDC instance.

Note that the open-access MAF files used for visualization have been processed with stringent filters to remove low quality and potential germline variants to ensure that only somatic variants are present. See (Masked Somatic Aggregation Workflow)[https://gdc-docs.nci.nih.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/#masked-somatic-aggregation-workflow] for details.

For additional reading and tutorials on the main features of cBioPortal please visit the main cBioPortal site hosted by [Memorial Sloan Kettering Cancer Center (MSKCC)](http://www.cbioportal.org/).

## Release Beta

* __GDC Product__: GDC cBioPortal
* __Release Date__: November 16, 2016


### Known Issues and Workarounds
*  Oncoprint feature is not active if multiple pipelines are selected <!--SV-512-->
*  Mutations classified as "Silent" or "RNA" in MAF file are not displayed in cBioPortal <!--SV-516-->
*  On study summary page the number of genes may not be equal on the gene list and bar chart.  This is because the gene list filters out non-cancer related genes with a single mutation in the selected project <!--SV-516-->
*  When using the browser's back button to return to the cBioPortal homepage, the "Select Patient/Case Set" menu may appear empty. Following the "To build your own case set" link generates an error message. Refreshing the page resolves the issue. <!-- SV-425 -->
