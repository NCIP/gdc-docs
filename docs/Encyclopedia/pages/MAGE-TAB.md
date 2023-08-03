# MAGE-TAB #
## Description ##

MicroArray Gene Expression Tabular (MAGE-TAB) archives are groups of tab-delimited spreadsheets that provide details about an experiment.  

## Overview ##

MAGE-TAB archives provide information about the data collection and data processing that was not performed by the GDC. The files in these archives can include molecular protocols, computational pipelines, or information about file organization. MAGE-TAB archives are available for download at the GDC Legacy Archive<sup>1</sup> with data files that were processed with some TCGA pipelines and can be downloaded from file pages or retrieved from the API (`.../files?expand=metadata_files`). The files in a MAGE-TAB archive will include information about all files in a complete study, not just for the associated data file.

### Structure ###

Contents of a MAGE-TAB file:

* Investigative Description Format (IDF): Provides general information about the study. This includes a brief description, the investigator's contact details, bibliographic references, and a text description of the protocols used in the study<sup>2</sup>
* Sample and Data Relationship Form (SDRF): Describes the relationships between samples, arrays, data, and other objects used or produced in the study<sup>2</sup>
* Array Design Format (ADF): Defines each array type used. An ADF file describes the design of an array, e.g., which sequence is located at each position on an array and associated annotations<sup>2</sup>
* Description File: Provides details about how the data files and molecular material were processed
* Changes File: Provides details about any changes that were made to each of the files in the MAGE-TAB archive
* Manifest File: Lists file names and md5sums of the files that should be included in the MAGE-TAB archive
* Readme File: Provides basic details about the MAGE-TAB archive and the associated study

## References ##
1. [GDC Legacy Archive](https://portal.gdc.cancer.gov/legacy-archive/)
2. [TCGA Encyclopedia - MAGE-TAB](https://wiki.nci.nih.gov/display/TCGA/MAGE-TAB)

## External Links ##
* N/A

Categories: Data Format
