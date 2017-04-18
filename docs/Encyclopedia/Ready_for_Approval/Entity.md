# Entity #
## Description ##
An entity in the GDC is a unique component of the GDC Data Model.

## Overview ##
The GDC Data Model is the primary method of organizing all data within the GDC<sup>1</sup>.  More specifically, all the data within the GDC can be thought of as a Directed Acyclic Graph (DAG) composed of interconnected entities.  A graphical representation of the GDC Data Model can be found [here](https://gdc.cancer.gov/developers/gdc-data-model/gdc-data-model-components).

Each entity in the GDC has a set of properties.  An example entity would be a `case` (patient). A `case` is linked to a number of other entities in the data model, such as those that contain Biospecimen and Clinical data.  An example entity that contains the ethnicity (*cases.demographic.ethnicity*) of the case is the `demographic` entity.  The GDC Data Model defines how each of the entities are connected and the GDC Data Dictionary defines the entities and the relationships between them<sup>2</sup>

Each entity is assigned a unique identifier in the form  of a version 4 UUID.

Data submitters can create and update submittable entities in the GDC Data Model and upload data files registered in the model using the GDC Data Submission Portal, the GDC API, and the GDC Data Transfer Tool<sup>3</sup>.

## References ##
1. [GDC Data Model](https://gdc.cancer.gov/developers/gdc-data-model)
2. [GDC Data Dictionary](https://docs.gdc.cancer.gov/Data_Dictionary/)
3. [GDC Data Submission Portal](https://gdc.cancer.gov/submit-data/gdc-data-submission-portal)

## External Links ##
* N/A
