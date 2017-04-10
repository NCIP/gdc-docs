# Biospecimen Data #

## Description ##
Biospecimen data describe biospecimen-related elements. A biospecimen element is a classification of biospecimen data that describe a common form of biospecimen. A biospecimen often takes the form of a case, sample, portion, analyte, or aliquot. They can be represented in a hierarchical relationship to one another.

## Overview ##

In the GDC, biospecimen elements refer to biological specimens that originate from cancer patients and each other.  For example, a `sample` element (entity) comes directly from a patient, but a `portion` that was used to generate a diagnostic slide comes from the `sample`. Each entity has associated attributes that can be used to describe it.  For example, an `analyte` entity has the `analyte_type` attribute, which reports the type of biological molecule that was extracted to produce the `analyte`.

Data included in the GDC will be associated with different biospecimen entities depending on their type. For example, clinical data files will be directly associated with a case entity (a patient) because they are associated with the patient as a whole. The data used to produce an RNA-Seq alignment originates from sequenced RNA extracts and would be associated with an aliquot entity. See the GDC Data Model and Data Dictionary for details on the biospecimen entity types and their properties.

The TCGA program, from which many of the GDC attributes were inherited, includes biospecimen elements that are reflected in the structure of the TCGA Barcode.

### Data ###

Biospecimen data can be downloaded from the GDC in several formats:

__API Retrieval:__ Specific information about each biospecimen entity can be queried from the API in a tab-delimited or JSON format. This information is available to retrieve programatically.
__Biospecimen Supplements:__ Biospecimen supplements are stored in XML format and can be downloaded from the GDC Portal. These XML files may contain biospecimen fields that do not appear in the API. These can be fields that only apply to certain projects or have not yet been incorporated into the GDC API.
__Biotab Files:__ Biotab files are available in the GDC Legacy Archive as tab-delimited files on a project-level basis. These may also include fields that are not available in the GDC API.

Additionally, the biospecimen data that is available from the API is displayed on the GDC Portal in the summary page for each case.

## References ##
1. [TCGA Encyclopedia - Biospecimen Data](https://wiki.nci.nih.gov/display/TCGA/Biospecimen+data)
2. [TCGA Encyclopedia - Biospecimen Element Type](https://wiki.nci.nih.gov/display/TCGA/Biospecimen+element+type)
3. [GDC Data Model](https://gdc.cancer.gov/developers/gdc-data-model/gdc-data-model-components)
4. [GDC Biospecimen Harmonization](https://gdc.cancer.gov/about-data/data-harmonization-and-generation/biospecimen-data-harmonization)
5. [GDC Data Dictionary](https://docs.gdc.cancer.gov/Data_Dictionary/viewer/)
6. [TCGA Barcode](https://wiki.nci.nih.gov/display/TCGA/TCGA+barcode)

## External Links ##
* N/A

Categories: Data Category
