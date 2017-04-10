# Clinical Supplement #

## Description ##

The purpose of a clinical supplement file is to make available clinical information that does not fit into an established data structure.

## Overview ##

The GDC Data Model has an extensive set of clinical elements, but the GDC model of data harmonization prohibits too many project-specific exceptions to be incorporated. This makes the Clinical Supplement files very useful as they may contain additional elements and are connected only to the associated case.    

For example, whether or not a patient has mutations in the ER, PR and HER2 genes is very important when researching breast cancer.  This information would not be very useful in most other cancers, so it can be found in the clinical supplement files for breast cancer patients rather than being incorporated into the data model across all cancers.

While some elements are unique to certain cancers, there is no guarantee that all patients associated with these cancers will have all elements populated. 

### Data Formats ###

Clinical Supplement files are available in the GDC Data Portal as case-level XML files. Most elements have an associated CDE-ID that can be queried at the CDE Browser website for additional information.  

## References ##
1. [GDC Clinical Data Elements](https://gdc.cancer.gov/about-data/data-harmonization-and-generation/clinical-data-harmonization)
2. [CDE Browser](https://cdebrowser.nci.nih.gov/cdebrowserClient/cdeBrowser.html#/search)

## External Links ##
* N/A

Categories: Data Type
