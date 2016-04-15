# Submission

## Overview

The GDC Submission API uses methods and endpoints that are distinct from those that drive the functionality of the GDC Data Portal. In particular, data and metadata that are in the process of being submitted can only be queried using [GraphQL](#querying-submitted-data-and-metadata-using-graphql), and not the methods described in [Search and Retrieval](Search_and_Retrieval.md).

This section describes the GDC API's submission functionality, including methods for submitting, deleting, updating, searching, and retrieving data and metadata.

## Submission endpoint

### Constructing the endpoint URL

The endpoint for submitting data to a specific project in the GDC is constructed as follows:
<pre>https://gdc-api.nci.nih.gov/<b>[&#x3C;API_version&#x3E;/]</b>submission/<b>&#x3C;Program.name&#x3E;</b>/<b>&#x3C;Project.code&#x3E;</b></pre>
where `[<API_version>/]` is the optional API version component (see [Getting Started](Getting_Started.md)).

#### Program.name and Project.code

The values of `Program.name` and `Project.code` can be obtained from the project URL on the GDC Data Submission Portal:

<pre>https://gdc-portal.nci.nih.gov/submission/<b>&#x3C;Program.name&#x3E;</b>/<b>&#x3C;Project.code&#x3E;</b>/dashboard</pre>

#### Dry Run Mode

The `submission` endpoint provides a `_dry_run` mode that simulates submission transactions without making changes to the GDC. This mode is activated by appending `/_dry_run` to the end of a submission endpoint.

#### Example

For example, a project with GDC Data Submission Portal URL

<pre>https://gdc-portal.nci.nih.gov/submission/<b>TCGA</b>/<b>ALCH</b>/dashboard</pre>

would have a versioned submission endpoint at

<pre>https://gdc-api.nci.nih.gov/<b>v0/</b>submission/<b>TCGA</b>/<b>ALCH</b></pre>

an unversioned submission endpoint at

<pre>https://gdc-api.nci.nih.gov/submission/<b>TCGA</b>/<b>ALCH</b></pre>

and a dry run endpoint at

<pre>https://gdc-api.nci.nih.gov/v0/submission/TCGA/ALCH/<b>_dry_run</b></pre>

[//]: # (this is just a comment ignore me I beg of you_)



## GDC Data Model

### Nodes, Properties, Links, and Schemas

The GDC Data Model is a representation of data stored in the GDC. It is used to retrieve, submit, update, and delete data. Although the GDC Data Model may contain some cyclic elements, it can be helpful to think of it as a [Directed Acyclic Graph (DAG)](https://en.wikipedia.org/wiki/Directed_acyclic_graph) composed of interconnected **nodes**. Each node in the GDC has a set of properties and links.

* **Properties** are key-value pairs associated with a node. Properties cannot be nested, which means that the value must be numerical, boolean, or a string, and cannot be another key-value set. Properties can be either required or optional. The following properties are of particular importance in constructing the GDC Data Model:
    * **Type** is a required property for all nodes. Node types include `project`, `case`, `demographic`, `sample`, `read_group` and others.
    * **System properties** are properties used in GDC system operation and maintenance, that cannot be modified except under special circumstances.
    * **Unique keys** are properties, or combinations of properties, that can be used to uniquely identify the node in the GDC. For example, the tuple (combination) of `[ project_id, submitter_id ]` is a unique key for most nodes, which means that although `submitter_id` does not need to be unique in GDC, it must be unique within a project.
* **Links** define relationships between nodes, and the multiplicity of those relationships (e.g. one-to-one, one-to-many, many-to-many).

The properties and links that a node can have are defined by the **JSON schema** corresponding to the node's `type`. Node JSON schemas are stored in the [GDC Data Dictionary](#gdc-data-dictionary).

Functionally similar node types are grouped under the same **category**. For example, node types `slide_image` and `submitted_unaligned_reads` belong to `data_file` category, which comprises nodes that correspond to files downloadable from the GDC Object Store.


### Unique Keys

When a node is created, it must be assigned a unique identifier in the form of a [version 4 universally unique identifier (UUID)](https://en.wikipedia.org/wiki/Universally_unique_identifier). The UUID uniquely identifies the node in the GDC, and is stored in the node's `id` property. For most submittable nodes, the UUID can be assigned by the submitter. If the submitter does not provide a UUID, it will be assigned by the GDC and returned in the API response upon successful completion of the transaction. See [Appendix D](Appendix_D_Format_of_Submission_Requests_and_Responses.md) for details of the API response format.

In addition to `id`, many nodes also include a `submitter_id` field. This field can contain any string (e.g. a "barcode") that the submitter wishes to use to identify the node. Typically this string identifies a corresponding entry in submitter's records. The GDC's only requirement with respect to `submitter_id` is that it be a string that is unique for all nodes within a project. The GDC Submission API requires a `submitter_id` for most nodes.

**Note:** For `case` nodes, `submitter_id` must correspond to a `submitted_subject_id` of a study participant registered with the project in dbGaP.

### GDC Data Dictionary

The [GDC Data Dictionary Viewer](../../Data_Dictionary/index.md) provides a user-friendly overview of node schemas, grouped by category. It also provides templates in JSON and TSV formats, that can be used for creating or updating nodes. Information in the GDC Data Dictionary can also be accessed programmatically using `dictionary` and `template` endpoints.

#### Dictionary Endpoint

The entire collection of GDC node schemas can be downloaded at the following GDC Data Dictionary endpoint:

<pre>https://gdc-api.nci.nih.gov/v0/submission/_dictionary/<b>_all</b></pre>

[//]: # (this is just a comment ignore me I beg of you_)

Individual schemas can be downloaded at the endpoint that corresponds to the node type. For example, the JSON schema for `case` nodes can be found at:

<pre>https://gdc-api.nci.nih.gov/v0/submission/_dictionary/<b>case</b></pre>

#### Template Endpoint

Submission templates are accessible programmatically at the `templates` endpoint. For example, the JSON template for `case` nodes can be obtained from:

	https://gdc-api.nci.nih.gov/v0/submission/template/case?format=json

and the entire collection of GDC Submission Templates can be obtained from:

	https://gdc-api.nci.nih.gov/v0/submission/template/?format=json


## Format of Submission API Requests and Responses

When creating or updating nodes in the GDC, the request must specify the node `type`, the node `id` or `submitter_id`, relationships (links) that the node has to existing nodes, and node properties as defined by the node JSON schema. To delete nodes, only the `id` property is required. The general format of GDC API submission requests and responses is provided in [Appendix D](Appendix_D_Format_of_Submission_Requests_and_Responses.md).

## Submission Transactions

Submission involves a series of transactions initiated by the submitter, that create and link nodes according to their [schemas](#nodes-properties-links-and-schemas), constructing a graph that can be represented by the diagram provided [here](https://gdc.nci.nih.gov/node/8396/). With the exception of `program` and `project` administrative nodes created by the GDC, all new nodes must be linked, at creation, to existing nodes or to new nodes being created in the same transaction. For example, a submitter cannot create a `portion` node unless the submitter either (1) has previously created the corresponding `case` and `sample` nodes, or (2) is creating those nodes in the same transaction. This also means that nodes cannot be deleted if they have "child" nodes attached to them.

If multiple nodes are being created and/or updated in a transaction, and an error is encountered for one of the nodes, then the transaction will fail and no changes will be made to the GDC.

## Creating and Updating Nodes

### POST and PUT Requests

The GDC Submission API provides two methods for creating nodes: HTTP POST requests and HTTP PUT requests:

* The **POST** method will create nodes that do not exist, and will fail if any of the nodes in the transaction already exist in the GDC.

* The **PUT** method will create new nodes and update existing nodes, and identify which nodes were created or updated in the API response.

The GDC suggests using POST for creating new nodes, and using PUT only for updating nodes. This helps to avoid inadvertent node updates that can occur when using PUT for creating nodes.

### Example: Creating and Updating a Case Node

In this example, a case node is created using POST. Then an attempt is made to create the same node again using POST, resulting in an error. Then the originally created node is updated (with the same information) using PUT.

The JSON in the request was generated using the `case` JSON template that can be obtained from the [GDC Data Dictionary Viewer](../../Data_Dictionary/index.md) and from `https://gdc-api.nci.nih.gov/v0/submission/template/case?format=json`.

**Note:** For `case` nodes, `submitter_id` must correspond to a `submitted_subject_id` of a study participant registered with the project in dbGaP.


```Request1
{
  "project_id": "TCGA-ALCH",
  "type": "case",
  "submitter_id": "TCGA-ALCH-000001",
  "projects": {
    "code": "ALCH"
  }

}
```
```Command1
export token=ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTO

curl --header "X-Auth-Token: $token" --request POST --data @Request https://gdc-api.nci.nih.gov/v0/submission/TCGA/ALCH
```
```Response1
{
  "cases_related_to_created_entities_count": 0,
  "cases_related_to_updated_entities_count": 0,
  "code": 201,
  "created_entity_count": 1,
  "entities": [
    {
      "action": "create",
      "errors": [],
      "id": "fbf69646-5904-4f95-92d6-692bde658f05",
      "related_cases": [],
      "type": "case",
      "unique_keys": [
        {
          "project_id": "TCGA-ALCH",
          "submitter_id": "TCGA-ALCH-000001"
        }
      ],
      "valid": true,
      "warnings": []
    }
  ],
  "entity_error_count": 0,
  "message": "Transaction successful.",
  "success": true,
  "transaction_id": 215,
  "transactional_error_count": 0,
  "transactional_errors": [],
  "updated_entity_count": 0
}
```
```Command2
curl --header "X-Auth-Token: $token" --request POST --data @Request https://gdc-api.nci.nih.gov/v0/submission/TCGA/ALCH
```
```Response2
{
  "cases_related_to_created_entities_count": 0,
  "cases_related_to_updated_entities_count": 0,
  "code": 400,
  "created_entity_count": 0,
  "entities": [
    {
      "action": null,
      "errors": [
        {
          "keys": [
            "id"
          ],
          "message": "Cannot create entity that already exists. Try updating entity (PUT instead of POST)",
          "type": "NOT_UNIQUE"
        }
      ],
      "id": null,
      "related_cases": [],
      "type": "case",
      "unique_keys": [
        {
          "project_id": "TCGA-ALCH",
          "submitter_id": "TCGA-ALCH-000001"
        }
      ],
      "valid": false,
      "warnings": []
    }
  ],
  "entity_error_count": 1,
  "message": "Transaction aborted due to 1 invalid entity.",
  "success": false,
  "transaction_id": null,
  "transactional_error_count": 0,
  "transactional_errors": [],
  "updated_entity_count": 0
}
```
```Command3
curl --header "X-Auth-Token: $token" --request PUT --data @Request https://gdc-api.nci.nih.gov/v0/submission/TCGA/ALCH
```
```Response3
{
  "cases_related_to_created_entities_count": 0,
  "cases_related_to_updated_entities_count": 0,
  "code": 200,
  "created_entity_count": 0,
  "entities": [
    {
      "action": "update",
      "errors": [],
      "id": "fbf69646-5904-4f95-92d6-692bde658f05",
      "related_cases": [],
      "type": "case",
      "unique_keys": [
        {
          "project_id": "TCGA-ALCH",
          "submitter_id": "TCGA-ALCH-000001"
        }
      ],
      "valid": true,
      "warnings": []
    }
  ],
  "entity_error_count": 0,
  "message": "Transaction successful.",
  "success": true,
  "transaction_id": 216,
  "transactional_error_count": 0,
  "transactional_errors": [],
  "updated_entity_count": 1
}
```



### Example: Creating an Aliquot Node

In this example, an `aliquot` node and a `sample` node are created in a single transaction. The `aliquot` is linked to `sample` which is linked to `case`. The first request is an example of using `submitter_id` properties to link nodes together. The second request is an example of using UUIDs for creating the links.

#### Request 1: Creating Links Using submitter_id

```Request
[
  {
    "type": "sample",
    "submitter_id": "TCGA-ALCH-000001-SAMPLE000001",
    "project_id": "TCGA-ALCH",
    "sample_type": "Primary Tumor",
    "sample_type_id": "01",
    "cases": {
      "submitter_id": "TCGA-ALCH-000001"
    }
  },
  {
    "type": "aliquot",
    "project_id": "TCGA-ALCH",
    "submitter_id": "TCGA-ALCH-000001-SAMPLE000001-ALIQUOT000001",
    "samples": {
      "submitter_id": "TCGA-ALCH-000001-SAMPLE000001"
    }
  }
]
```
```Command
export token=ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTO

curl --header "X-Auth-Token: $token" --request POST --data @Request https://gdc-api.nci.nih.gov/v0/submission/TCGA/ALCH
```
```Response
{
  "cases_related_to_created_entities_count": 1,
  "cases_related_to_updated_entities_count": 0,
  "code": 201,
  "created_entity_count": 2,
  "entities": [
    {
      "action": "create",
      "errors": [],
      "id": "48270338-6464-448f-bbef-b09d4f80b11b",
      "related_cases": [
        {
          "id": "fbf69646-5904-4f95-92d6-692bde658f05",
          "submitter_id": "TCGA-ALCH-000001"
        }
      ],
      "type": "sample",
      "unique_keys": [
        {
          "project_id": "TCGA-ALCH",
          "submitter_id": "TCGA-ALCH-000001-SAMPLE000001"
        }
      ],
      "valid": true,
      "warnings": []
    },
    {
      "action": "create",
      "errors": [],
      "id": "7af58da0-cb3e-43e2-a074-4bd8f27565ba",
      "related_cases": [
        {
          "id": "fbf69646-5904-4f95-92d6-692bde658f05",
          "submitter_id": "TCGA-ALCH-000001"
        }
      ],
      "type": "aliquot",
      "unique_keys": [
        {
          "project_id": "TCGA-ALCH",
          "submitter_id": "TCGA-ALCH-000001-SAMPLE000001-ALIQUOT000001"
        }
      ],
      "valid": true,
      "warnings": []
    }
  ],
  "entity_error_count": 0,
  "message": "Transaction successful.",
  "success": true,
  "transaction_id": 222,
  "transactional_error_count": 0,
  "transactional_errors": [],
  "updated_entity_count": 0
}
```

#### Request 2: Creating Links Using UUID


```Request
[
  {
    "type": "sample",
    "submitter_id": "TCGA-ALCH-000001-SAMPLE000001",
    "id": "2aa7a07b-e706-4eef-aeba-b849972423a0",
    "project_id": "TCGA-ALCH",
    "sample_type": "Primary Tumor",
    "sample_type_id": "01",
    "cases": {
      "id": "fbf69646-5904-4f95-92d6-692bde658f05"
    }
  },
  {
    "type": "aliquot",
    "project_id": "TCGA-ALCH",
    "submitter_id": "TCGA-ALCH-000001-SAMPLE000001-ALIQUOT000001",
    "samples": {
      "id": "2aa7a07b-e706-4eef-aeba-b849972423a0"
    }
  }
]
```
```Command
export token=ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTO

curl --header "X-Auth-Token: $token" --request POST --data @Request https://gdc-api.nci.nih.gov/v0/submission/TCGA/ALCH
```
```Response
{
  "cases_related_to_created_entities_count": 1,
  "cases_related_to_updated_entities_count": 0,
  "code": 201,
  "created_entity_count": 2,
  "entities": [
    {
      "action": "create",
      "errors": [],
      "id": "2aa7a07b-e706-4eef-aeba-b849972423a0",
      "related_cases": [
        {
          "id": "fbf69646-5904-4f95-92d6-692bde658f05",
          "submitter_id": "TCGA-ALCH-000001"
        }
      ],
      "type": "sample",
      "unique_keys": [
        {
          "project_id": "TCGA-ALCH",
          "submitter_id": "TCGA-ALCH-000001-SAMPLE000001"
        }
      ],
      "valid": true,
      "warnings": []
    },
    {
      "action": "create",
      "errors": [],
      "id": "545096d5-ce1c-433f-80f0-fd0b04b56cb6",
      "related_cases": [
        {
          "id": "fbf69646-5904-4f95-92d6-692bde658f05",
          "submitter_id": "TCGA-ALCH-000001"
        }
      ],
      "type": "aliquot",
      "unique_keys": [
        {
          "project_id": "TCGA-ALCH",
          "submitter_id": "TCGA-ALCH-000001-SAMPLE000001-ALIQUOT000001"
        }
      ],
      "valid": true,
      "warnings": []
    }
  ],
  "entity_error_count": 0,
  "message": "Transaction successful.",
  "success": true,
  "transaction_id": 219,
  "transactional_error_count": 0,
  "transactional_errors": [],
  "updated_entity_count": 0
}
```


## Retrieving Nodes

JSON objects representing submitted nodes can be retrieved using the `entities` endpoint of the GDC Submission API. This endpoint retrieves nodes by UUID. See [GraphQL](#graphql) for more advanced methods of querying submitted data.

```Command
export token=ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTO

curl --header "X-Auth-Token: $token" https://gdc-api.nci.nih.gov/v0/submission/TCGA/ALCH/entities/fbf69646-5904-4f95-92d6-692bde658f05
```
```Response
{
  "entities": [
    {
      "program": "TCGA",
      "project": "ALCH",
      "properties": {
        "created_datetime": "2016-04-14T08:44:43.361800-05:00",
        "id": "fbf69646-5904-4f95-92d6-692bde658f05",
        "project_id": "TCGA-ALCH",
        "projects": [
          {
            "id": "d9906779-f1da-5d9f-9caa-6d5ecb2e3cd6",
            "submitter_id": null
          }
        ],
        "state": "validated",
        "submitter_id": "TCGA-ALCH-000001",
        "type": "case",
        "updated_datetime": "2016-04-14T21:29:28.401212-05:00"
      }
    }
  ]
}
```

## Deleting Nodes

The `entities` endpoint is also be used to delete nodes. This is accomplished using a DELETE request to the endpoint, specifying the node's UUID. If a node cannot be deleted because it is linked to child nodes, the GDC Submission API will respond with an error providing a list of nodes that must be deleted prior to deleting the subject node.


```Shell
export token=ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTO

curl --header "X-Auth-Token: $token" --request DELETE https://gdc-api.nci.nih.gov/v0/submission/TCGA/ALCH/entities/67782964-0065-491d-b051-2ae404bb734d
```
```Response
{
  "code": 200,
  "deleted_entity_count": 1,
  "dependent_ids": "",
  "entities": [
    {
      "action": "delete",
      "errors": [],
      "id": "67782964-0065-491d-b051-2ae404bb734d",
      "related_cases": [],
      "type": "case",
      "valid": true,
      "warnings": []
    }
  ],
  "entity_error_count": 0,
  "message": "Successfully deleted 1 entities",
  "success": true,
  "transaction_id": 192,
  "transactional_error_count": 0,
  "transactional_errors": []
}
```


## Uploading Data Files

Experimental data files like BAM and FASTQ can be uploaded directly to the API using the `files` endpoint, by specifying the UUID of the corresponding `data_file` node:

	export token=ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTO

	curl --request PUT --header "X-Auth-Token: $token" https://gdc-api.nci.nih.gov/v0/submission/TCGA/ALCH/files/6d45f2a0-8161-42e3-97e6-e058ac18f3f3 -d@data.fastq



## Querying Submitted Data Using GraphQL

**NOTE:** The GDC Submission API GraphQL service is an authenticated
  resource for which a GDC Authorization Token must be
  provided. Access is limited to authorized submitters.

### GraphQL Overview

"GraphQL is a query language designed to build client applications by
providing an intuitive and flexible syntax and system for describing
their data requirements and interactions." (from [GraphQL specification](https://facebook.github.io/graphql/)).

GraphQL has proven to be a very effective method of querying the hierarchical
nature of the GDC's graph datamodel.  The `/graphql` endpoint on the
GDC Submission API provides a real-time view of the state of the
entities in a project.

### Sample GraphQL query

The following is a GraphQL query for a case in project
TCGA-LAML that returns a JSON document containing the `submitter_id` of
the case and of its samples.


```GraphQL
{
    case (project_id: "TCGA-LAML", first: 1) {
         submitter_id
         samples { submitter_id }
    }
}
```
```Response
{
  "data": {
    "case": [
      {
        "samples": [
          {
            "submitter_id": "TCGA-AB-2901-11A"
          },
          {
            "submitter_id": "TCGA-AB-2901-03A"
          }
        ],
        "submitter_id": "TCGA-AB-2901"
      }
    ]
  }
}
```

### GDC Data Dictionary Usage

All fields defined in the GDC Data Dictionary can be queried using GraphQL.  The
GraphQL schema is generated off of the Data Dictionary.  For example,
if the term `submitter_id` was changed to `alias` for all Sample
Entities, the above query would be updated to contain `samples { alias
}` rather than `samples { submitter_id }`.

### GraphiQL IDE

GDC provides an [in-browser IDE](https://gdc-portal.nci.nih.gov/submission/graphiql) for exploring
GraphQL. It is an instance of [GraphiQL](https://github.com/graphql/graphiql).

This IDE provides tab-completion and syntax checking using the GraphQL
schema generated from the GDC Data Dictionary.  GraphiQL allows for
easy discoverability of both fields and query filters.


### API Usage

All authorized submitters (those who have READ permissions on a
project) can access project Entities via the GraphQL endpoint, located
at `/submission/graphql/`.  As with all authorized API requests, the
authorization token (as downloaded from the Portal) is provided in the
request header `X-Auth-Token`

The following example demonstrates a download using `curl` in a Unix
environment with one assumption that users have their token stored in the
environment variable `TOKEN`:

```bash
$ curl -XPOST -H"X-Auth-Token: $TOKEN" "https://gdc-api.nci.nih.gov/v0/submission/graphql/" -d'{"query": "{ aliquot(first: 2) { id }}"}'
```
```Response
{
  "data": {
    "aliquot": [
      {
        "id": "222417ce-fd47-4dcc-b458-e06751097099"
      },
      {
        "id": "7646a246-5b14-4a50-9c12-7c10dba580d1"
      }
    ]
  }
}
```

**NOTE:** Query results will only contain results from
  the projects that the user has READ access to.

**NOTE:** Query results have a default limit of 10 results, to choose
  a different number of results, override the default with `first: X`
  where `X` is the maximum number of desired results.  If `X` is `0`,
  then no limit is applied. (For pagination, see the `offset` argument)


### Examples

#### Example

Complete list of existing cases in a project titled "GDC-EXAMPLE", including count of cases and fields `submitter_id` and `id`:

```Query
{
  case (project_id: "GDC-EXAMPLE", first: 0) {
    id
    submitter_id

  }
  _case_count (project_id: "GDC-EXAMPLE")
}
```
```Response
{
  "data": {
    "_case_count": 20,
    "case": [
      {
        "id": "700d1110-b6b4-4251-89d4-fa6f0698e3f8",
        "submitter_id": "GDC-EXAMPLE-000004"
      },
      {
        "id": "be01357d-7348-40b4-a997-8a61ae7af17d",
        "submitter_id": "GDC-EXAMPLE-000005"
      },
      {
        "id": "e5638697-6ef3-4bf8-a373-102519093f33",
        "submitter_id": "GDC-EXAMPLE-000008"
      },
      {
        "id": "4871d41a-680e-4fd0-901c-b06f06ecae33",
        "submitter_id": "GDC-EXAMPLE-000007"
      },
      {
        "id": "2f18c2c1-bff2-43b6-9702-e138c72d8c6b",
        "submitter_id": "GDC-EXAMPLE-000009"
      },
      {
        "id": "ec83e038-4f01-47a6-bc69-47fb297d0282",
        "submitter_id": "GDC-EXAMPLE-000006"
      },
      {
        "id": "e4642952-d259-4be1-9c53-ed95aa1fc50b",
        "submitter_id": "GDC-EXAMPLE-000011"
      },
      {
        "id": "8bcaf0b3-21d0-45c6-87ee-c997efb417dc",
        "submitter_id": "GDC-EXAMPLE-000010"
      },
      {
        "id": "83de027e-bcbf-4239-975b-7e8ced82448e",
        "submitter_id": "GDC-EXAMPLE-000013"
      },
      {
        "id": "bbd91cc1-06e2-4e60-8b93-e09c3b16f00c",
        "submitter_id": "GDC-EXAMPLE-000014"
      },
      {
        "id": "574fd163-4368-440c-9548-d76a0fbc9056",
        "submitter_id": "GDC-EXAMPLE-000015"
      },
      {
        "id": "47c92cdd-ff11-4c25-b0f0-0f7671144271",
        "submitter_id": "GDC-EXAMPLE-000016"
      },
      {
        "id": "9f13caab-1fda-4b2a-b500-f79dc978c6c1",
        "submitter_id": "GDC-EXAMPLE-000017"
      },
      {
        "id": "9418f194-8741-44db-bd8f-36f4fd8c3bf2",
        "submitter_id": "GDC-EXAMPLE-000018"
      },
      {
        "id": "6fb2a018-c5f3-45e5-81d3-e58e7e4bf921",
        "submitter_id": "GDC-EXAMPLE-000019"
      },
      {
        "id": "70236972-e796-414a-9b7a-3b29b849ba7c",
        "submitter_id": "GDC-EXAMPLE-000020"
      },
      {
        "id": "6f78e86f-9e31-4af5-a0d9-b8970ece476d",
        "submitter_id": "GDC-EXAMPLE-000021"
      },
      {
        "id": "c6fcb2f0-c6bb-4b40-a761-bae3e63869cb",
        "submitter_id": "GDC-EXAMPLE-000002"
      },
      {
        "id": "67782964-0065-491d-b051-2ae404bb734d",
        "submitter_id": "GDC-EXAMPLE-000001"
      },
      {
        "id": "b45d2891-ba81-4ecc-a250-c58060934227",
        "submitter_id": "GDC-EXAMPLE-000012"
      }
    ]
  }
}
```
#### Example

GraphQL query for any one case in 'TCGA-LUAD' without Diagnosis information

```JavaScript
{
  case (project_id: "TCGA-LUAD", without_links: ["diagnoses"], first: 1) {
      submitter_id
  }
}
```
```Response
{
  "data": {
    "case": [
      {
        "submitter_id": "TCGA-17-Z050"
      }
    ]
  }
}
```

#### Example

GraphQL query for the number of cases in 'TCGA-LUAD' without Diagnosis information

```JavaScript
{
  _case_count (project_id: "TCGA-LUAD", without_links: ["diagnoses"])
}
```
```Response
{
  "data": {
    "_case_count": 5
  }
}
```

#### Example

Query for the release state of aliquots belonging to case with `submitter_id:
"TCGA-17-Z050"`

```JavaScript
{
  aliquot(with_path_to: {type: "case", submitter_id:"TCGA-17-Z050"}) {
    id release_state
  }
}
```

#### Example

GraphQL query that uses a graphql fragment to get specific properties from two portions
and give them aliases in the response.

```JavaScript
{
  some_portion: portion (first: 1) {
    ...portionProperties
  }
  specific_portion: portion(submitter_id: "TCGA-67-6217-01A-13-2191-20") {
    ...portionProperties
  }
}

fragment portionProperties on portion {
  submitter_id
  is_ffpe
}
```
```Response
{
  "data": {
    "some_portion": [
      {
        "is_ffpe": false,
        "submitter_id": "TCGA-62-A471-10A-01"
      }
    ],
    "specific_portion": [
      {
        "is_ffpe": false,
        "submitter_id": "TCGA-67-6217-01A-13-2191-20"
      }
    ]
  }
}
```

#### Example

GraphQL Query for a case in "TCGA-LUAD" and return a biospecimen tree

```JavaScript
{
  case(project_id: "TCGA-LUAD", first: 1) {
    id
    samples(first: 1) {
      id
      portions(first: 1) {
        id
        analytes(first: 1) {
          id
          aliquots(first: 1) {
            id
          }
        }
      }
    }
  }
}
```
```Response
{
  "data": {
    "case": [
      {
        "id": "19ca36e6-2154-4224-89b1-117a4a4407f6",
        "samples": [
          {
            "id": "5e2625d2-290d-48cd-af5c-27dc8e3c8b6a",
            "portions": [
              {
                "analytes": [
                  {
                    "aliquots": [
                      {
                        "id": "8e1820d5-dcd8-4760-9962-221e2b71d4b9"
                      }
                    ],
                    "id": "6449533c-e52a-4e58-bae7-0732f48153ef"
                  }
                ],
                "id": "26b75643-8fcd-445e-a0e0-9868cac589ea"
              }
            ]
          }
        ]
      }
    ]
  }
}
```
