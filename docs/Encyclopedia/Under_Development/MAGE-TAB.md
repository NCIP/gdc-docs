# MAGE-TAB #
## Description ##

MicroArray Gene Expression Tabular (MAGE-TAB) archives are groups of tab-delimited spreadsheets that provide details about an experiment.  

## Overview ##

MAGE-TAB archives are available for download at the GDC Legacy Archive with data files that were processed with some TCGA pipelines. The information in a MAGE-TAB archive will include information about all files in a complete study, not just for the associated data file.

### Structure ###

Contents of a MAGE-TAB file:

* Investigative Description Format: Provides general information about the investigation. This includes its name, a brief description, the investigator's contact details, bibliographic references, and free text descriptions of the protocols used in the investigation<sup>XXX</sup>.
* Sample and Data Relationship Form: Describes the relationships between samples, arrays, data, and other objects used or produced in the investigation. Each row represents an analyzed element (often an aliquot) in its most basic electronic form (i.e. raw data file) and the production of higher-level data files as protocols (e.g. normalization) are applied to the file and its derivatives. These protocols correspond to those listed in the IDF.
* Array Design Format: Defines each array type used. An ADF file describes the design of an array, e.g., what sequence is located at each position on an array and what the annotation of this sequence is.
* Description File: Provides details about how the data files and molecular material were processed
* Changes File: Provides details about any changes that were made to the files in the MAGE-TAB archive
* Manifest File: Lists file names and md5sums of the files that should be included in the MAGE-TAB archive
* Readme File: Provides basic details about the MAGE-TAB archive and the associated study.

## References ##
1. [TCGA Encyclopedia - MAGE-TAB](https://wiki.nci.nih.gov/display/TCGA/MAGE-TAB)

## External Links ##
* TBD

Categories: Data Format
