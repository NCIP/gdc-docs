# GDC Data Dictionary

## Introduction

The GDC Data Dictionary is a file that describes the content of a database’s metadata, data about data, and how the data will be used by a system. It defines the structure of a database, the [data model](../Data/Data_Model/GDC_Data_Model.md), and the rules in which the data in a system needs to abide by. The data dictionary contains records about other objects in the database, such as data ownership and data relationships to other objects.

### A Data Dictionary Includes the Following:

* Comprehensive list of names of tables in the database and their schemas.
* Comprehensive list of indexes and columns to which the tables in those indexes relate.
* Constraints defined on tables, including primary keys, foreign-key relationships to other tables, and not-null constraints.

### Standards and Conventions

The Center for Translational Data Science (CTDS) is expanding the information found in the GDC Data Dictionary by including references to external standards such as the National Cancer Institute Thesaurus (NCIt) and the Cancer Data Standards Registry and Repository (caDSR). These databases share overlapping domains with terms related to cancer, drugs, and temporal events. The ability to cross reference against multiple databases will provide researchers with properly documented references and compatible terminologies across studies. This enforcement of congruency within the documentation of metadata allows for a more robust harmonization of data within the GDC.

## Data Dictionary Viewer

The [GDC Data Dictionary Viewer](viewer.md) is a user-friendly interface for accessing the dictionary. It includes the following functionality:

*   __Dictionary contents:__ Display of entities defined in the dictionary, including their descriptions, properties, and links.
*   __Links to semantic resources:__ Links to semantic data resources that define [Common Data Elements (CDEs)](http://cde.nih.gov) used in the dictionary
*   __Submission templates:__ JSON and TSV template generation for use in GDC data submission.

### Components of the Data Dictionary Viewer

#### Title and Summary

[![Title and Summary](images/GDC_DD_Title_and_Summary.png)](images/GDC_DD_Title_and_Summary.png "Click to see the full image.")

* __Type:__ The name of the node.
* __Category:__ The type of metadata; some examples are Clinical, Biospecimen, Analysis and Submittable Data Files.
* __Description:__ This section contains a written explanation for the type of data that would be found in this node.
* __Unique Keys:__ The properties that are specific to this node.

This section also contains a "Download Template" link with a drop down menu containing the two template file types, TSV and JSON. These files will contain all properties that are found in the node, but [not all properties are required](#properties) to upload the node. 

#### Links

[![Links](images/GDC_DD_Links.png)](images/GDC_DD_Links.png "Click to see the full image.")

* __Links to Entity:__ Other nodes that are connected to the active node.
* __Link Name:__ The name of the node as it would be used when referring to one of its fields in the active node. Example: Refering to the case's `submitter_id` field in the sample node, the link would be `cases.submitter_id`.
* __Relationship:__ The written description for the association between the active node and the other connected node.
* __Required:__ Displays whether the node is required for the existence of the active node.

#### Properties

[![Properties Enumeration](images/GDC_DD_Properties_Enumeration.png)](images/GDC_DD_Properties_Enumeration.png "Click to see the full image.")
[![Properties Integer](images/GDC_DD_Properties_Integer.png)](images/GDC_DD_Properties_Integer.png "Click to see the full image.")
[![Properties Number](images/GDC_DD_Properties_Number.png)](images/GDC_DD_Properties_Number.png "Click to see the full image.")
[![Properties String](images/GDC_DD_Properties_String.png)](images/GDC_DD_Properties_String.png "Click to see the full image.")
[![Properties Boolean](images/GDC_DD_Properties_Boolean.png)](images/GDC_DD_Properties_Boolean.png "Click to see the full image.")

* __Property:__ The field name found in the node.

* __Description:__ The written explanation for the expected type and characterization of data found in this field.

* __Acceptable Values:__ The values that are expected to be entered into the field based on the value category:
    * Enumeration: A list of predetermined strings. The user must select the exact string, case matters, from the list to be a valid entry. Many of these properties with enumerations have numerous values, to see all of the values, click the "More Values" link at the bottom of the property row under the __Acceptable Values__ column.
    * Integer: A field that only accepts whole numbers.
    * Number: A field that can accept any number, even numbers with decimal places. 
    * String: A field in which alphanumeric characters and `_`, `.`, `-`, up to a length of 32,767, can be entered. Do not use other characters as it will create submission errors.
    * Boolean: A field that only accpets `true` or `false` as acceptable values. For a boolean value, the `true` and `false` values are all lowercase. If these values are not entered as lowercase, the dictionary will not recognize the value and create an error.

* __Required:__ This informs the user whether this field is necessary for the submission of the node.

* __CDE:__ The CDE Public ID, with the direct link to its respective Data Element Details page.


## GDC Metadata Validation Service

The MVS tool enables easier query of the GDC Data Dictionary for data submitters and recommends GDC properties and values based on user-supplied synonyms.  Created by the NCI CBIIT EVS Team, it leverages NCI vocabulary systems caDSR and NCI Thesaurus. Below are some of the features included in the MVS tool:

*   Users can complete partial or exact match searches
*   Searches can include terms that are synonymous to the GDC allowable values
*   Users can compare their list of values to the GDC allowable values
*   Dictionary paths are described so users can find the specific node where a property is located

## Entity JSON Schemas

In technical terms, the dictionary is a set of YAML files that define JSON schemas for each entity in the dictionary. The files are available [on GitHub](https://github.com/NCI-GDC/gdcdictionary/tree/develop/gdcdictionary/schemas).

The GDC API can generate entity JSON schemas in JSON format. The API also provides the template generation functionality accessible via the GDC Data Dictionary Viewer. See [API documentation](../API/Users_Guide/Submission/#gdc-data-dictionary-endpoints) for details on how to access these functions programmatically.
