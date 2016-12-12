# Mutation Annotation Format (MAF) #
## Introduction ##
A Mutation Annotation Format (MAF) is a tab-delimited file containing somatic and/or germline mutation annotations. MAF files containing any germline mutation annotations are kept in the controlled access portion of the GDC Data Portal. MAF files containing only somatic mutations are kept in the open access portion of the Data Portal.[1]
## Description ##
### Overview ###
Mutations are discovered by aligning DNA sequences derived from tumor samples to sequences derived from normal samples and a reference sequence.  A MAF file identifies, for each sample, the discovered putative or validated mutations and categorizes those mutations (SNP, deletion, or insertion) as somatic (originating in the tissue) or germline (originating from the germline) as well as the annotation for those mutations.[1]

A MAF file identifies, for each sample, the discovered putative or validated mutations and categorizes those mutations (SNP, deletion, or insertion) as somatic (originating in the tissue) or germline (originating from the germline). These can be subcategorized as follows:

*Somatic mutations:*
  - Missense and nonsense
  - Splice site, defined as SNP within 2 bp of the splice junction
  - Silent mutations
  - Indels that overlap the coding region or splice site of a gene or the targeted region of a genetic element of interest.
  - Frameshift mutations
  - Mutations in regulatory regions
  
*SNPs:*
  - Any germline SNP with validation status "unknown" is included.
  - SNPs already validated in dbSNP are not included since they are unlikely to be involved in cancer.[1]

### Structure ###
The structure of the MAF is available in the [TCGA MAF Specification] (https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification). Changes made to the MAF in support of the GDC are available in the [GDC MAF Format document] (https://gdc-docs.nci.nih.gov/Data/File_Formats/MAF_Format/).
## Resources ##
| Resource | Location |
| TCGA MAF specification | https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification |
| GDC MAF Format | https://gdc-docs.nci.nih.gov/Data/File_Formats/MAF_Format/ | 
## References ##
[1] [Mutation Annotation Format] (https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format)
