# Biospecimen Data #

## Description ##
Biospecimen data includes information about how a patient's tissue was processed and subsampled for use in histology and molecular assays. They can be represented in a hierarchical relationship to one another.

## Overview ##

In the GDC, biospecimen elements refer to biological specimens that originate from cancer patients. For example, a `sample` element (entity) comes directly from a patient and a `portion` that was used to generate a diagnostic slide comes from the `sample`. Each entity has associated attributes that can be used to describe it. For example, an `analyte` entity has the `analyte_type` attribute, which reports the type of biological molecule that was extracted to produce the `analyte`.

Data files hosted at the GDC will be associated with different biospecimen entities depending on their type. For example, clinical data files will be directly associated with a case entity (a patient) because these values are associated with the patient as a whole. The data used to produce an RNA-Seq alignment originates from sequenced RNA extracts and would be associated with an aliquot entity. See the GDC Data Model<sup>1</sup> and Data Dictionary<sup>2</sup> for details on the biospecimen entity types and their properties.

The TCGA program, from which many of the GDC attributes were inherited, includes biospecimen elements that are reflected in the structure of the TCGA Barcode<sup>3</sup>.

### Data ###

Biospecimen data can be downloaded from the GDC in several formats:

* __API Retrieval:__ Specific information about each biospecimen entity can be queried from the API in a tab-delimited or JSON format. This information can be retrieved programatically.
* __Biospecimen Supplements:__ Biospecimen supplements are stored in XML format and can be downloaded from the GDC Portal. Supplements may contain biospecimen fields that do not appear in the API as well as those that do. These can include fields that only apply to certain projects or have not yet been incorporated into the GDC Data Dictionary.
* __Biotab Files:__ Biotab files specific supplemental files that are available in the GDC Legacy Archive as tab-delimited files on a project-level basis. These may also include fields that are not available in the GDC API.

Additionally, the biospecimen data that is available from the API is displayed on the GDC Portal in the summary page for each case.

## References ##
1. [GDC Data Dictionary](https://docs.gdc.cancer.gov/Data_Dictionary/viewer/)
2. [GDC Data Model](https://gdc.cancer.gov/developers/gdc-data-model/gdc-data-model-components)
3. [TCGA Barcode](https://wiki.nci.nih.gov/display/TCGA/TCGA+barcode)

## External Links ##
* [GDC Biospecimen Harmonization](https://gdc.cancer.gov/about-data/data-harmonization-and-generation/biospecimen-data-harmonization)

Categories: Data Category
