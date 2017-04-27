# MAGE-TAB #
## Description ##

MicroArray Gene Expression Tabular (MAGE-TAB) archives are groups of tab-delimited spreadsheets that provide details about an experiment.  

## Overview ##

MAGE-TAB archives provide information about data collection and data processing that was not performed by the GDC.

MAGE-TAB archives are available for download at the GDC Legacy Archive<sup>1</sup> with data files that were processed with some TCGA pipelines. The files in a MAGE-TAB archive will include information about all files in a complete study, not just for the associated data file.

### Structure ###

Contents of a MAGE-TAB file:

* __Investigative Description Format (IDF):__ Provides general information about the investigation. This includes its name, a brief description, the investigator's contact details, bibliographic references, and free text descriptions of the protocols used in the investigation<sup>2</sup>.
* __Sample and Data Relationship Form (SDRF):__ Describes the relationships between samples, arrays, data, and other objects used or produced in the investigation. Each row represents an analyzed element (often an aliquot) in its most basic electronic form (i.e. raw data file) and the production of higher-level data files as protocols (e.g. normalization) are applied to the file and its derivatives. These protocols correspond to those listed in the IDF<sup>2</sup>.
* __Array Design Format (ADF):__ Defines each array type used. An ADF file describes the design of an array, e.g., which sequence is located at each position on an array and associated annotations<sup>2</sup>.
* __Description File:__ Provides details about how the data files and molecular material were processed
* __Changes File:__ Provides details about any changes that were made to the files in the MAGE-TAB archive
* __Manifest File:__ Lists file names and md5sums of the files that should be included in the MAGE-TAB archive
* __Readme File:__ Provides basic details about the MAGE-TAB archive and the associated study

## References ##
1. [GDC Legacy Archive](https://portal.gdc.cancer.gov/legacy-archive/)
2. [TCGA Encyclopedia - MAGE-TAB](https://wiki.nci.nih.gov/display/TCGA/MAGE-TAB)

## External Links ##
* N/A

Categories: Data Format
