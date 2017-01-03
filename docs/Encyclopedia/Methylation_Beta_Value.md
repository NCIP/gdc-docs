# Methylation Beta Value #
## Description ##
A Methylation Beta Value represents the ratio between the methylated and total array intensities.  
## Overview ##

The Methylation Beta Value files at the GDC associate methylation beta values with chromosomal regions and biological features of that region.   


### Data Formats ###

Methylation Beta Value files are available as tab-delimited plain text files.  See the GDC Methylation Liftover documentation for descriptions of each field.  

```
Composite Element REF	Beta_value	Chromosome	Start	End	Gene_Symbol	Gene_Type	Transcript_ID	Position_to_TSS	CGI_Coordinate	Feature_Type
cg00000292	0.484837056351681	chr16	28878779	28878780	ATP2A1;ATP2A1;ATP2A1;ATP2A1;ATP2A1	protein_coding;protein_coding;protein_coding;protein_coding;protein_coding	ENST00000357084.6;ENST00000395503.7;ENST00000536376.4;ENST00000562185.4;ENST00000563975.1	373;290;-1275;-465;-83	CGI:chr16:28879633-28880547	N_Shore
cg00002426	0.314288770620393	chr3	57757816	57757817	SLMAP;SLMAP;SLMAP;SLMAP;SLMAP;SLMAP	protein_coding;protein_coding;protein_coding;protein_coding;protein_coding;protein_coding	ENST00000295951.6;ENST00000295952.6;ENST00000383718.6;ENST00000428312.4;ENST00000449503.5;ENST00000467901.1	1585;368;261;257;257;514	CGI:chr3:57756198-57757263	S_Shore
cg00003994	0.607727818648732	chr7	15686237	15686238	MEOX2	protein_coding	ENST00000262041.5	576	CGI:chr7:16399497-16399700	.
```

## References ##
1. [GDC Methylation Liftover Pipeline Documentation](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Methylation_LO_Pipeline/)
2. Zhou, Wanding, Laird Peter L., and Hui Shen. "Comprehensive characterization, annotation and innovative use of Infinium DNA methylation BeadChip probes." Nucleic Acids Research. (2016): doi: 10.1093/nar/gkw967

## External Links ##
* TBD

Categories: Data Type
