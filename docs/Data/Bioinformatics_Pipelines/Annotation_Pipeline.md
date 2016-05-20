# GDC Somatic Mutation Annotation Pipeline

GDC annotates raw VCF files using [Variant Effect Predictor (VEP)](http://www.ensembl.org/info/docs/tools/vep/index.html) v84 ([reference](http://dx.doi.org/10.1093/bioinformatics/btq330)) and GDC plugins, including the following additional datasets:

*   COSMIC v.75
*   GENCODE v.22
*   sift v.5.2.2
*   ESP v.20141103
*   polyphen v.2.2.2
*   dbSNP v.146
*   Ensembl genebuild v.2014-07
*   Ensembl regbuild v.13.0
*   HGMD public v.20154
*   ClinVar v.201601


See [GDC VCF Format: GDC INFO fields](../File_Formats/VCF_Format.md#gdc-info-fields) for information on variant annotation fields that are currently included in Annotated Somatic Mutation VCF files.
