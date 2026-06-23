# Case #

## Description ##
A case is the collection of all data related to a specific patient in the context of a specific project.

## Overview ##

A case refers to a specific cancer patient or cell line in the context of a project<sup>1</sup>. In the GDC Data Model, cases are associated with a project, the samples collected from the case, and the clinical data<sup>2</sup>.

Cases are associated with all clinical data elements as well as the clinical and biospecimen supplemental information<sup>3</sup>. The case entity also serves as an API endpoint which allows for case-level data to be easily retrievable using programmatic methods.

Cases are also the basic unit of registration for the GDC Submission Portal. All cases that are to be submitted to a project need to be registered with dbGaP on a case-level basis prior to submission.  

## References ##
1. [GDC Data Dictionary - Case](https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#?view=table-definition-view&id=case)
2. [GDC Data Model](https://gdc.cancer.gov/developers/gdc-data-model/gdc-data-model-components)
3. [GDC Portal Case Detail Page](https://docs.gdc.cancer.gov/Data_Portal/Users_Guide/Repository/#cases-list)


Categories: General
