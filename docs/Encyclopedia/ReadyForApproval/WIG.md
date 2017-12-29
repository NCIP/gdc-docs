# Wiggle (WIG) #

## Description ##
A text data format containing alignment scores over genomic intervals.

## Overview ##

The WIG (wiggle) format is designed for display of dense continuous data such as probability scores. Wiggle data elements must be equally sized; if you need to display continuous data that is sparse or contains elements of varying size, use the BedGraph format instead.

After alignment, sequence reads are typically summarized into scores over/within genomic intervals.  The resulting files are formatted differently:
* BED files - genomic intervals with additional information
* Wiggle files, BEDgraphs, BigWigs - genomic intervals with scores
* GFF/GTF - genomic annotation with information and scores

The wiggle file format has two formats: variableStep and fixedStep.

### variableStep Wiggle Format ###

variableStep format is designed for data with irregular intervals between data points, and is the more commonly used format. It begins with a declaration line, followed by two columns containing chromosome positions and data values.

The declaration line begins with the word variableStep and is followed by space-separated key-value pairs:

* __chrom (required)__ - name of chromosome
* __span (optional, defaults to 1)__ - the number of bases that each data value should cover
The span allows data to be compressed as follows:

Without span:

variableStep chrom=chr2
300701  12.5
300702  12.5
300703  12.5
300704  12.5
300705  12.5

With span:

variableStep chrom=chr2 span=5
300701  12.5
Both of these examples will display a value of 12.5 at position 300701-300705 on chromosome 2.

### fixedStep Wiggle Format ###

fixedStep format is designed for data with regular intervals between data points, and is the more compact of the two wiggle formats. It begins with a declaration line, followed by a single column of data values.

The declaration line begins with the word fixedStep and is followed by space-separated key-value pairs:

* __chrom (required)__ - name of chromosome
* __start (required)__ - start point for the data values
* __step (required)__ - distance between data values
* __span (optional, defaults to 1)__ - the number of bases that each data value should cover
Without span

fixedStep chrom=chr3 start=400601 step=100
11
22
33

Displays the values 11, 22, 33 as single-base features, on chromosome 3 at positions 400601, 400701 and 400801 respectively.

With span

fixedStep chrom=chr3 start=400601 step=100 span=5
11
22
33

## References ##
1. N/A

## External Links ##
* [Ensembl Wiggle Format Definition](http://useast.ensembl.org/info/website/upload/wig.html)

Categories: Data Format
