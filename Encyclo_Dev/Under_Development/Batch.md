# Batch #

## Description ##

A batch describes a group of [aliquots](Aliquot.md) from which data was collected in the same experiment using the same protocols.  More specifically, a batch can refer to a group of aliquots from which data was produced in the same sequencing flow cell or on the same array chip.  

Batch effects are artifactual similarities between datasets that were produced simultaneously. Batch effects come from the method of data collection rather than the genetics of the sample in question.

## Overview ##

In the GDC, batches can be determined in several ways. Array-based data files that come from TCGA start with a five letter code that corresponds to an array plate. Additionally, TCGA barcodes have an embedded code that corresponds to the plate in which the aliquot was processed.   

## References ##
1. [TCGA Barcode](https://wiki.nci.nih.gov/display/TCGA/TCGA+barcode)

## External Links ##
* N/A

Categories: General
