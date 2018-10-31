# GDC Data Dictionary

## Introduction

The GDC Data Dictionary defines components of the [GDC Data Model](../Data/Data_Model/GDC_Data_Model.md) and relationships between them.

## Data Dictionary Viewer

The [GDC Data Dictionary Viewer](viewer.md) is a user-friendly interface for accessing the dictionary. It includes the following functionality:

*   _Dictionary contents:_ Display of entities defined in the dictionary, including their descriptions, properties, and links.
*   _Links to semantic resources:_ Links to semantic data resources that define [Common Data Elements (CDEs)](http://cde.nih.gov) used in the dictionary
*   _Submission templates:_ Generation JSON and TSV templates for use in GDC data submission.

## GDC Metadata Validation Service

The MVS tool enables easier query of the GDC Data Dictionary for data submitters and recommends GDC properties and values based on user-supplied synonyms.  Created by the NCI CBIIT EVS Team, it leverages NCI vocabulary systems caDSR and NCI Thesaurus. Below are some of the features included in the MVS tool:

*   Users can complete partial or exact match searches
*   Searches can include terms that are synonymous to the GDC allowable values
*   Users can compare their list of values to the GDC allowable values
*   Dictionary paths are described so users can find the specific node where a property is located

## Entity JSON Schemas

In technical terms, the dictionary is a set of YAML files that define JSON schemas for each entity in the dictionary. The files are available [on GitHub](https://github.com/NCI-GDC/gdcdictionary/tree/develop/gdcdictionary/schemas).

The GDC API can generate entity JSON schemas in JSON format. The API also provides the template generation functionality accessible via the GDC Data Dictionary Viewer. See [API documentation](../API/Users_Guide/Submission/#gdc-data-dictionary-endpoints) for details on how to access these functions programmatically.
