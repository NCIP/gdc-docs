# Copy Number Segment #

## Description ##

Copy-number segment files are generated from the GDC [Copy Number Variation](LINK) Pipeline using the [DNAcopy](LINK) R-package.  

## Overview ##

The main function of the harmonized Copy Number Segment files at the GDC is to identify contiguous segments of the chromosome that had undergone duplications or deletions. This is performed by locating probes with adjacent binding regions with intensities that suggest a
similar copy number.    


### Data Formats ###

Copy Number Segment files are tab-delimited and associate chromosomal coordinates with copy number data including log<sub>2</sub> ratio segment means and the number of probes associated with that segment.  

```
Sample	Chromosome	Start	End	Num_Probes	Segment_Mean
6_A01_466074	1	61735	668210	33	0.1974
6_A01_466074	1	690090	7406118	3568	-0.0755
6_A01_466074	1	7407908	7512202	82	-0.8006
6_A01_466074	1	7517702	12020317	2412	-0.0949
```



## References ##
1.[GDC CNV Analysis Pipeline Documentation](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/CNV_Pipeline/)

## External Links ##
* TBD

Categories: Data Type
