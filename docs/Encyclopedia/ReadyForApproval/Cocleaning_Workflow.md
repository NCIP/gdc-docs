# Cocleaning Workflow #

## Description ##

Cocleaning is the process in which Whole Exome Sequencing alignments are optimized by comparing paired tumor/normal samples. The WXS BAM files that are downloadable from the GDC Portal have been co-cleaned after alignment [1].  

## Overview ##

The DNA-Seq alignment cocleaning workflow comprises two steps:
* Base quality score recalibration (BQSR):  Base quality scores are reassessed based on errors across both samples.
* Indel realignment: Insertion-deletion mutations are realigned to minimize the number of base mismatches.

### Tools ###

Both cocleaning steps are implemented using two functions of the Genome Analysis Toolkit (GATK): BaseRecalibrator and IndelRealigner [2].  

## References ##
1. [GDC DNA-Seq Pipeline Documentation](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/)
2. [GATK Website](https://software.broadinstitute.org/gatk/)

## External Links ##
* N/A

Categories: Workflow Type
