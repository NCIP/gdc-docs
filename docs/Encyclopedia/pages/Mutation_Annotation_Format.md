# Mutation Annotation Format (MAF) #
## Description ##
Mutation Annotation Format (MAF) files are tab-delimited files that contain somatic and/or germline mutation annotations. MAF files containing any germline mutation annotations are protected and distributed in the controlled access portion of the GDC Data Portal. MAF files containing only somatic mutations are publicly available.

## Overview ##
Variants are discovered by aligning DNA sequences derived from tumor samples and sequences derived from normal samples to a reference sequence<sup>1</sup>. A MAF file identifies, for all samples in a project, the discovered putative or validated mutations and categorizes those mutations (polymorphism, deletion, or insertion) as somatic (originating in the tumor tissue) or germline (originating from the germline). Mutation annotation are also reported. Note that VCF files report on all transcripts affected by a mutation while MAF files only report on the most affected one.  

MAF files are generated at the GDC using the Variant Aggregation pipeline<sup>2</sup> by aggregating case-level VCF files on a project and pipeline level (one protected MAF per project per pipeline). Protected MAF files are further processed into open-access somatic MAF files by removing low quality or germline mutations.

### Structure ###
The structure of the MAF is available in the [GDC MAF Specification](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/) along with descriptions for each field<sup>3</sup>.

## References ##
1. [GDC DNA-Seq Analysis](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/)
2. [GDC-Dictionary: Somatic Aggregation Pipeline](https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=somatic_aggregation_workflow)
3. [GDC MAF TCGA (Legacy)Format ](https://docs.gdc.cancer.gov/Encyclopedia/pages/Mutation_Annotation_Format_TCGA/)



Categories: Data Format
