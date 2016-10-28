# Submission

## Overview

The GDC Submission API uses methods and endpoints that are distinct from those that drive the functionality of the GDC Data Portal. In particular, data and metadata that are in the process of being submitted can only be queried using [GraphQL](#querying-submitted-data-and-metadata-using-graphql), and not the methods described in [Search and Retrieval](Search_and_Retrieval.md).

This section describes the GDC API's submission functionality, including methods for submitting, deleting, updating, searching, and retrieving data and metadata.

## Submission endpoint

### Constructing the endpoint URL

The endpoint for submitting data to a specific project in the GDC is constructed as follows:
<pre>https://gdc-api.nci.nih.gov/<b>[&#x3C;API_version&#x3E;/]</b>submission/<b>&#x3C;Program.name&#x3E;</b>/<b>&#x3C;Project.code&#x3E;</b></pre>
where `[<API_version>/]` is the optional API version component (see [Getting Started](Getting_Started.md)).

The values of `Program.name` and `Project.code` can be obtained from the project URL on the GDC Data Submission Portal:

<pre>https://gdc-portal.nci.nih.gov/submission/<b>&#x3C;Program.name&#x3E;</b>/<b>&#x3C;Project.code&#x3E;</b>/dashboard</pre>

For more information about program name and project code see [The GDC Data Model  section](../../Data/Data_Model/GDC_Data_Model/#program-name-project-code-and-project-id).

#### Example

For example, a project with GDC Data Submission Portal URL

<pre>https://gdc-portal.nci.nih.gov/submission/<b>TCGA</b>/<b>ALCH</b>/dashboard</pre>

would have a versioned submission endpoint at

<pre>https://gdc-api.nci.nih.gov/<b>v0/</b>submission/<b>TCGA</b>/<b>ALCH</b></pre>

and an unversioned submission endpoint at

<pre>https://gdc-api.nci.nih.gov/submission/<b>TCGA</b>/<b>ALCH</b></pre>


## GDC Data Model

Submitters should review the [GDC Data Model documentation](../../Data/Data_Model/GDC_Data_Model.md) and the [GDC Data Dictionary](../../Data_Dictionary/index.md) before initiating submission.

### UUIDs

Submitters can assign UUIDs to all submittable entities other than those that correspond to user-downloadable files. If the submitter does not provide a UUID, it will be assigned by the GDC and returned in the API response upon successful completion of the submission transaction. See [Appendix C](Appendix_C_Format_of_Submission_Requests_and_Responses.md) for details of the API response format. To learn more about UUIDs see the [GDC Data Model documentation](../../Data/Data_Model/GDC_Data_Model.md#uuids).

### Submitter IDs

In addition to `id`, many entities also include a `submitter_id` field. This field can contain any string (e.g. a "barcode") that the submitter wishes to use to identify the entity. Typically this string identifies a corresponding entry in submitter's records. The GDC's only requirement with respect to `submitter_id` is that it be a string that is unique for all entities within a project. The GDC Submission API requires a `submitter_id` for most entities.

**Note:** For `case` entities, `submitter_id` must correspond to a `submitted_subject_id` of a study participant registered with the project in dbGaP.

### GDC Data Dictionary Endpoints

Information in the [GDC Data Dictionary](../../Data_Dictionary/index.md) can be accessed programmatically as described below.

#### Submission Templates

Submission templates are accessible programmatically at the `templates` endpoint. Template format (`json`, `tsv` or `csv`) is specified using the `format` parameter.

For example, the JSON template for `case` entities can be obtained from:

	https://gdc-api.nci.nih.gov/v0/submission/template/case?format=json

A set of templates for all entities in the GDC Data Model can be downloaded from:

	https://gdc-api.nci.nih.gov/v0/submission/template/?format=json

#### Entity JSON Schemas

The entire collection of GDC entity schemas can be downloaded from the `dictionary` endpoint:

<pre>https://gdc-api.nci.nih.gov/v0/submission/_dictionary/<b>_all</b></pre>

[//]: # (this is just a comment ignore me I beg of you_)

Individual schemas can be downloaded by specifying entity type. For example, the JSON schema for `case` entities can be found at:

<pre>https://gdc-api.nci.nih.gov/v0/submission/_dictionary/<b>case</b></pre>


## Format of Submission API Requests and Responses

When creating or updating entities in the GDC, the request must specify the entity `type`, the entity `id` or `submitter_id`, relationships (links) that the entity has to existing entities, and entity properties as defined by the [GDC Data Dictionary](../../Data_Dictionary/index.md). To delete entities, only the `id` property is required. The general format of GDC API submission requests and responses is provided in [Appendix C](Appendix_C_Format_of_Submission_Requests_and_Responses.md).

## Submission Transactions

Submission of data to the GDC involves a series of transactions initiated by the submitter, that create and link entities according to the [GDC Data Model](../../Data/Data_Model/GDC_Data_Model.md). With the exception of `program`, which is an administrative entity created by the GDC, all new entities must be linked, at creation, to existing entities or to new entities being created in the same transaction. For example, a submitter cannot create a `portion` entity unless the submitter either (1) has previously created the corresponding `case` and `sample` entities, or (2) is creating those entities in the same transaction. This also means that entities cannot be deleted if they have "child" entities attached to them.

If multiple entities are being created and/or updated in a transaction, and an error is encountered for one of the entities, then the transaction will fail and no changes will be made to the GDC.

### Dry Run Transactions

The `submission` endpoint provides a `_dry_run` mode that simulates submission transactions without making changes to the GDC. This mode is activated by appending `/_dry_run` to the end of a submission endpoint.

The following is an example of a POST request, that simulates creating an entity in dry run mode:


```Request
{
  "project_id": "TCGA-ALCH",
  "type": "case",
  "submitter_id": "TCGA-ALCH-000001",
  "projects": {
    "code": "ALCH"
  }
}
```
```Command
export token=ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTO

curl --header "X-Auth-Token: $token" --request POST --data @Request https://gdc-api.nci.nih.gov/v0/submission/TCGA/ALCH/_dry_run
```
```Response
{
  "cases_related_to_created_entities_count": 0,
  "cases_related_to_updated_entities_count": 0,
  "code": 200,
  "created_entity_count": 1,
  "entities": [
    {
      "action": "create",
      "errors": [],
      "id": "61f48d1c-9439-448c-a90c-d6dbe76b3654",
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
  "message": "Transaction would have been successful. User selected dry run option, transaction aborted, no data written to database.",
  "success": true,
  "transaction_id": null,
  "transactional_error_count": 0,
  "transactional_errors": [],
  "updated_entity_count": 0
}
```

#### Dry Run Commit

For convenience, the GDC enables users to commit earlier `_dry_run` transactions instead of uploading the same data again to execute the changes. This `commit` action is allowed on transactions that (1) have not been previously committed and (2) were successful `dry_run` transactions.

Note that the `commit` action is a separate transaction with its own transaction id, and it can be executed [asynchronously](#asynchronous-transactions). If the state of the submission project has changed in a way that would make the original `_dry_run` transaction invalid if it were run again (e.g. entities with the same `submitter_id` have since been created in another transaction), then then `commit` action will fail.

```Command
export token=ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTO

curl --header "X-Auth-Token: $token" --request POST https://gdc-api.nci.nih.gov/v0/submission/TCGA/ALCH/transactions/467/commit?async=true
```
```Response
{
  "code": 200,
  "message": "Transaction submitted.",
  "transaction_id": 468,
}
```



#### Dry Run Close

The GDC Submission API also provides a `close` action on `_dry_run` transactions. This `close` action is allowed on `_dry_run` transactions that have not been previously closed. Closing a `_dry_run` transaction prevents it from being committed in the future.

Note that the `close` action is a separate transaction with its own transaction id.

```Command
export token=ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTO

curl --header "X-Auth-Token: $token" --request POST https://gdc-api.nci.nih.gov/v0/submission/TCGA/ALCH/transactions/467/close
```
```Response
{
    "code": 200,
    "message": "Closed transaction.",
    "transaction_id": <transaction_id>
}
```




### Asynchronous Transactions

The `submission` endpoint provides an asynchronous mode that provides immediate response and executes submission transactions in the background. This mode is activated by appending `?async=true` to the end of a submission endpoint.  The API will respond with the `transaction_id` which can be used to look up the result of the transaction at a later time via the [GraphQL](#querying-submitted-data-and-metadata-using-graphql) endpoint.  If the server has too many asynchronous jobs scheduled already, your request to schedule a transaction may fail.

#### Example

The following is an example of a PUT request, that creates a case asynchronously:

```Request
{
  "project_id": "TCGA-ALCH",
  "type": "case",
  "submitter_id": "TCGA-ALCH-000001",
  "projects": {
    "code": "ALCH"
  }
}
```
```Command
export token=ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTO

curl --header "X-Auth-Token: $token" --request POST --data @Request https://gdc-api.nci.nih.gov/v0/submission/TCGA/ALCH?async=true
```
```Response
{
  "code": 200,
  "message": "Transaction submitted.",
  "transaction_id": 467,
}
```

The following is a [GraphQL](#querying-submitted-data-and-metadata-using-graphql) request that looks up the state of the above transaction:

```GraphQL_Request
query {
  transaction_log(id: 467) {
    is_dry_run
    committed_by
    state
  }
}
```
```GraphQL_Response
{
  "data": {
    "transaction_log": [
      {
        "committed_by": null,
        "is_dry_run": false,
        "state": "FAILED"
      }
    ]
  }
}
```

#### Transaction Status

The following transaction fields can be queried using [GraphQL](#querying-submitted-data-and-metadata-using-graphql) and are helpful in determining the status of a transaction:

|Field|Type|Description|
|---|---|---|
|`id`|ID|Transaction identifier|
|`is_dry_run`|Boolean|Indicates whether the transaction is a dry run|
|`closed`|Boolean|For dry run transactions, indicates whether the transaction has been closed to prevent it from being committed in the future.|
|`committable`|Boolean|Indicates whether the transaction can be committed (i.e. it is a successful dry run transaction that has not been committed previously and has not been closed)|
|`state`|String|Indicates the state of the transaction: `PENDING`, `SUCCEEDED`, `FAILED` (due to user error), or `ERRORED` (due to system error)|
|`committed_by`|ID|The ID of the transaction that committed this transaction|

**Note:** To check whether a dry run transaction was committed successfully, check the `state` of the transaction that executed the commit. The `state` of the dry run transaction itself does not represent the status of a subsequent commit.

## Creating and Updating Entities

### POST and PUT Requests

The GDC Submission API provides two methods for creating entities: HTTP POST requests and HTTP PUT requests:

* The **POST** method will create entities that do not exist, and will fail if any of the entities in the transaction already exist in the GDC.

* The **PUT** method will create new entities and update existing entities, and identify which entities were created or updated in the API response.

The GDC suggests using POST for creating new entities, and using PUT only for updating entities. This helps to avoid inadvertent entity updates that can occur when using PUT for creating entities.



### Example: Creating and Updating Case Entities

In this example, a case entity is created using POST. Then an attempt is made to create the same entity again using POST, resulting in an error. Then the originally created entity is updated (with the same information) using PUT.

The JSON in the request was generated using the `case` JSON template that can be obtained from the [GDC Data Dictionary Viewer](../../Data_Dictionary/index.md) and from `https://gdc-api.nci.nih.gov/v0/submission/template/case?format=json`.

**Note:** For `case` entities, `submitter_id` must correspond to a `submitted_subject_id` of a study participant registered with the project in dbGaP.


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



### Example: Creating an Aliquot Entity

In this example, an `aliquot` entity and a `sample` entity are created in a single transaction. The `aliquot` is linked to `sample` which is linked to `case`. The first request is an example of using `submitter_id` properties to link entities together. The second request is an example of using UUIDs for creating the links.

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


## Retrieving Entities

### Entities Endpoint

JSON objects representing submitted entities can be retrieved using the `entities` endpoint of the GDC Submission API. This endpoint retrieves entities by UUID. A single UUID or a comma-separated list of UUIDs can be passed to this endpoint as a query.

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

### Export Endpoint

The `export` endpoint provides additional functionality for exporting entities from the GDC submission system. The `ids` parameter accepts a UUID or a comma-separated list of UUIDs. The `format` parameter allows the user to specify the preferred format of the API response: JSON, TSV, or CSV. When the `with_children` parameter is set to `with_children`, the response includes the metadata stored in all "child" entities of the entity being requested.


```Command
export token=ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTO


curl --header "X-Auth-Token: $token" 'https://gdc-api.nci.nih.gov/v0/submission/TCGA/ALCH/export?ids=11f8321-832f-4a8b-8384-a2f6256557e0&format=json&with_children=with_children'
```
```Response
{
  "case": [
    {
      "tissue_source_sites": [],
      "submitter_id": "TCGA-ALCH-000026",
      "project_id": "TCGA-ALCH",
      "type": "case",
      "id": "11f83251-832f-4a8b-8384-a2f6256557e0",
      "projects": [
        {
          "code": "ALCH",
          "id": "d9906779-f1da-5d9f-9caa-6d5ecb2e3cd6"
        }
      ]
    }
  ],
  "sample": [
    {
      "sample_type_id": "10",
      "time_between_excision_and_freezing": null,
      "oct_embedded": "false",
      "tumor_code_id": null,
      "submitter_id": "Blood-00001_api26",
      "intermediate_dimension": null,
      "id": "23308708-6a63-471e-947c-6a93c6e85983",
      "time_between_clamping_and_freezing": null,
      "pathology_report_uuid": null,
      "tumor_descriptor": null,
      "sample_type": "Blood Derived Normal",
      "project_id": "TCGA-ALCH",
      "current_weight": null,
      "composition": null,
      "is_ffpe": null,
      "shortest_dimension": null,
      "tumor_code": null,
      "tissue_type": null,
      "days_to_sample_procurement": null,
      "cases": [
        {
          "id": "11f83251-832f-4a8b-8384-a2f6256557e0",
          "submitter_id": "TCGA-ALCH-000026"
        }
      ],
      "freezing_method": null,
      "type": "sample",
      "preservation_method": null,
      "days_to_collection": null,
      "initial_weight": null,
      "longest_dimension": null
    }
  ],
  "read_group": [
    {
      "library_name": "Solexa-34688",
      "is_paired_end": true,
      "size_selection_range": null,
      "adapter_sequence": null,
      "library_strand": null,
      "submitter_id": "Blood-00001-aliquot_lane1_barcode26",
      "library_preparation_kit_name": null,
      "adapter_name": null,
      "target_capture_kit_name": null,
      "includes_spike_ins": null,
      "library_preparation_kit_version": null,
      "id": "90163202-cfd7-4f6a-8214-e7e4e924d3a6",
      "spike_ins_concentration": null,
      "target_capture_kit_vendor": null,
      "read_length": 75,
      "sequencing_date": "2010-08-04",
      "spike_ins_fasta": null,
      "to_trim_adapter_sequence": null,
      "RIN": null,
      "platform": "Illumina",
      "library_selection": "Hybrid_Selection",
      "library_strategy": "WXS",
      "library_preparation_kit_catalog_number": null,
      "target_capture_kit_target_region": null,
      "fastq_name": null,
      "target_capture_kit_version": null,
      "aliquots": [
        {
          "id": "e66dee54-5f4c-4471-9e08-dba0f6cdaaa4",
          "submitter_id": "Blood-00001-aliquot26"
        }
      ],
      "read_group_name": "205DD.3-2",
      "library_preparation_kit_vendor": null,
      "project_id": "TCGA-ALCH",
      "type": "read_group",
      "target_capture_kit_catalog_number": null,
      "instrument_model": "Illumina HiSeq 2000",
      "base_caller_name": null,
      "experiment_name": "Resequencing",
      "flow_cell_barcode": "205DDABXX",
      "sequencing_center": "BI",
      "base_caller_version": null
    }
  ],
  "aliquot": [
    {
      "source_center": "23",
      "centers": [],
      "analytes": [],
      "submitter_id": "Blood-00001-aliquot26",
      "amount": 10,
      "samples": [
        {
          "id": "23308708-6a63-471e-947c-6a93c6e85983",
          "submitter_id": "Blood-00001_api26"
        }
      ],
      "concentration": 0.07,
      "project_id": "TCGA-ALCH",
      "type": "aliquot",
      "id": "e66dee54-5f4c-4471-9e08-dba0f6cdaaa4"
    }
  ],
  "submitted_unaligned_reads": [
    {
      "read_groups": [
        {
          "id": "90163202-cfd7-4f6a-8214-e7e4e924d3a6",
          "submitter_id": "Blood-00001-aliquot_lane1_barcode26"
        }
      ],
      "data_type": "Unaligned Reads",
      "file_name": "dummy.fastq",
      "md5sum": "70c48a8a670ed2a02327601a10038d06",
      "data_format": "FASTQ",
      "submitter_id": "Blood-00001-aliquot_lane1_barcode26.fastq",
      "state_comment": null,
      "data_category": "Sequencing Data",
      "file_size": 38,
      "project_id": "TCGA-ALCH",
      "type": "submitted_unaligned_reads",
      "id": "6d45f2a0-8161-42e3-97e6-e058ac18f3f3",
      "experimental_strategy": "WGS"
    },
    {
      "read_groups": [
        {
          "id": "90163202-cfd7-4f6a-8214-e7e4e924d3a6",
          "submitter_id": "Blood-00001-aliquot_lane1_barcode26"
        }
      ],
      "data_type": "Unaligned Reads",
      "file_name": "dummy.fastq",
      "md5sum": "70c48a8a670ed2a02327601a10038d06",
      "data_format": "FASTQ",
      "submitter_id": "Blood-00001-aliquot_lane1_barcode27.fastq",
      "state_comment": null,
      "data_category": "Sequencing Data",
      "file_size": 38,
      "project_id": "TCGA-ALCH",
      "type": "submitted_unaligned_reads",
      "id": "4faabdd6-45bb-4259-8868-13d5b1149748",
      "experimental_strategy": "WGS"
    }
  ]
}
```

### GraphQL

Submitters can use the GraphQL query language for advanced search and retrieval of data from the GDC Submission Portal. See [GraphQL](#querying-submitted-data-using-graphql) for more information.


## Deleting Entities

The `entities` endpoint can also be used to delete entities. This is accomplished using a DELETE request to the endpoint, specifying the entity's UUID. If an entity cannot be deleted because it is linked to child entities, the GDC Submission API will respond with an error providing a list of entities that must be deleted prior to deleting the subject entity.

A subgraph (a parent along with all of its child entities) can be deleted in a single transaction by passing a comma-separated list of UUIDs to the `entities` endpoint.

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

## Working With Files

### Uploading Data Files

Experimental data files like BAM and FASTQ can be uploaded directly to the API using the `files` endpoint, by specifying the UUID of the corresponding `data_file` entity.  Uploading files may be more efficiently performed using the [GDC Data Transfer Tool](/Data_Transfer_Tool/Users_Guide/Getting_Started.md).

	export token=ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTO

	curl --request PUT --header "X-Auth-Token: $token" https://gdc-api.nci.nih.gov/v0/submission/TCGA/ALCH/files/6d45f2a0-8161-42e3-97e6-e058ac18f3f3 -d@data.fastq

#### Upload Manifest

The `manifest` endpoint generates a manifest for uploading files using the GDC Data Transfer Tool. It requires a comma-separated list of file UUIDs to generate a manifest.

	https://gdc-api.nci.nih.gov/v0/submission/PROGRAM/PROJECT/manifest?ids=bf0751ca-fc3b-4760-b876-0fefce040be5,90163202-cfd7-4f6a-8214-e7e4e924d3a6

### Downloading Files

Unreleased files that have been uploaded to the GDC can be downloaded by submitters using the `data` endpoint and an appropriate authentication token. See [Downloading Files](Downloading_Files.md) for details.

### Deleting Files

Uploaded files can be deleted by deleting the entity that corresponds to the file. See [Deleting Entities](#deleting-entities) for details.

## Querying Submitted Data Using GraphQL

### GraphQL Overview

[GraphQL](https://facebook.github.io/graphql/) is a query language that makes it easy to search and retrieve data from graph data structures such as the GDC Data Model.

Unlike the methods outlined in [Search and Retrieval](Search_and_Retrieval.md), which provide access to public releases (or snapshots) of GDC data, the `/graphql` endpoint of GDC Submission API makes it possible for submitters to access "live" data, which provides a real-time view of the state of entities in a project.

**NOTE:** Access to GDC Submission API GraphQL service is limited to authorized and authenticated submitters. Submitters may only access data in their own project using GraphQL.


### GraphQL IDE

The GDC GraphQL IDE is an instance of [GraphiQL](https://github.com/graphql/graphiql), an in-browser GraphQL IDE that facilitates construction and execution of GraphQL queries. The GDC GraphQL IDE provides tab-completion and syntax checking using schema from the GDC Data Dictionary. It can be found at [https://gdc-portal.nci.nih.gov/submission/graphiql](https://gdc-portal.nci.nih.gov/submission/graphiql).

Before interacting directly with the GDC Submission API's GraphQL endpoint, users are encouraged to become familiar with executing queries using the GDC GraphQL IDE.

### GraphQL Endpoint

GDC data submitters can access the GDC Submission API GraphQL endpoint at:

<pre>https://gdc-api.nci.nih.gov/[&#x3C;API_version&#x3E;/]submission/&#x3C;Program.name&#x3E;/&#x3C;Project.code&#x3E;<b>/graphql</b></pre>

**NOTE:** An authentication token is required for all requests to the `graphql` endpoint. Queries are restricted to those projects for which the submitter has obtained authorization.


### Constructing a Query

When sending GraphQL requests to the API directly, the bare GraphQL query must be wrapped in a "query" JSON object as shown below:

<pre>
{
	"query": "<b>{Bare_GraphQL_Query}</b>",
	"variables": null
}
</pre>

When using the GDC GraphQL IDE, the bare JSON query must be used without a JSON wrapper.

#### Bare GraphQL query

In its simplest form, a GraphQL query is a **selection set** (curly brackets) that encloses a set of **fields**. The selection set defines the set of information that is to be retrieved. Furthermore, in GraphQL fields are conceptually equivalent to functions that retrieve additional fields and, in some cases, can take arguments. So each field in a selection set can have its own selection set, thereby creating a nested query structure that can navigate complex data relationships. See [GraphQL Specification](https://facebook.github.io/graphql/) for further details.

In GDC GraphQL IDE, a root field (field within the outermost/umbrella selection set) typically corresponds to an entity, whereas fields inside nested selection sets are typically a combination of entities and entity properties.

The "Docs" panel on the right-hand side of the GDC GraphQL IDE allows users to discover the fields that can be queried with GraphQL. Note that the panel contains a lot of information and users may experience a delay before it is displayed.

A simple GraphQL query looks like this:

	{
	  case (project_id: "TCGA-ALCH", first: 0) {
	    id
	    submitter_id

	  }
	  _case_count (project_id: "TCGA-ALCH")
	}

[//]: # (this is just a comment ignore me I beg of you_)


The query above has two root fields: `case` and `_case_count`. The `case` field corresponds to the `case` entity in the GDC Data Model. The query supplies two arguments to the field:

1. `project_id: "TCGA-ALCH"`, which requests only cases in the TCGA-ALCH project.
2. `first: 0`, which requests that the API provide all results in the response, without pagination ( a nonzero positive integer value of `first` specifies the number of results to return, 10 by default; "pages" are selected using `offset`).

The `_case_count` field is a special field that returns the number of cases that match the supplied argument.

The bare query above can be used as is in the GraphQL IDE. In order to pass this query to the GDC API directly, it needs to be further processed as described below.

#### Passing GraphQL queries to GDC API directly

Before a bare GraphQL query is passed to the GDC API, it must be processed as follows:

1. [Escape](http://text-rescue.com/string-escape/json-escape-tool.html) the query using JSON string rules
2. Wrap the query in a ["query" JSON object](#constructing-a-query).
3. Pass the query to the `graphql` endpoint in an HTTP POST request.

Using the `case` and `_case_count` example above as the starting point, the results are as follows:

```bare_GraphQL
{
	case (project_id: "TCGA-ALCH", first: 0) {
		id
		submitter_id

	}
	_case_count (project_id: "TCGA-ALCH")
}
```
```escaped_GraphQL
{\n\tcase (project_id: \"TCGA-ALCH\", first: 0) {\n\t\tid\n\t\tsubmitter_id\n\n\t}\n\t_case_count (project_id: \"TCGA-ALCH\")\n}
```
```Query_json
{
	"query": "{\n\tcase (project_id: \"TCGA-ALCH\", first: 0) {\n\t\tid\n\t\tsubmitter_id\n\n\t}\n\t_case_count (project_id: \"TCGA-ALCH\")\n}",
	"variables": null
}
```
```Shell_command
export token=ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTO

curl --request POST --header "X-Auth-Token: $token" 'https://gdc-api.nci.nih.gov/v0/submission/graphql' -d@Query_json
```
```API_Response
{
  "data": {
    "_case_count": 20,
    "case": [
      {
        "id": "700d1110-b6b4-4251-89d4-fa6f0698e3f8",
        "submitter_id": "TCGA-ALCH-000004"
      },
      {
        "id": "be01357d-7348-40b4-a997-8a61ae7af17d",
        "submitter_id": "TCGA-ALCH-000005"
      },
      {
        "id": "e5638697-6ef3-4bf8-a373-102519093f33",
        "submitter_id": "TCGA-ALCH-000008"
      },
      {
        "id": "4871d41a-680e-4fd0-901c-b06f06ecae33",
        "submitter_id": "TCGA-ALCH-000007"
      },
      {
        "id": "2f18c2c1-bff2-43b6-9702-e138c72d8c6b",
        "submitter_id": "TCGA-ALCH-000009"
      },
      {
        "id": "ec83e038-4f01-47a6-bc69-47fb297d0282",
        "submitter_id": "TCGA-ALCH-000006"
      },
      {
        "id": "e4642952-d259-4be1-9c53-ed95aa1fc50b",
        "submitter_id": "TCGA-ALCH-000011"
      },
      {
        "id": "8bcaf0b3-21d0-45c6-87ee-c997efb417dc",
        "submitter_id": "TCGA-ALCH-000010"
      },
      {
        "id": "83de027e-bcbf-4239-975b-7e8ced82448e",
        "submitter_id": "TCGA-ALCH-000013"
      },
      {
        "id": "bbd91cc1-06e2-4e60-8b93-e09c3b16f00c",
        "submitter_id": "TCGA-ALCH-000014"
      },
      {
        "id": "574fd163-4368-440c-9548-d76a0fbc9056",
        "submitter_id": "TCGA-ALCH-000015"
      },
      {
        "id": "47c92cdd-ff11-4c25-b0f0-0f7671144271",
        "submitter_id": "TCGA-ALCH-000016"
      },
      {
        "id": "9f13caab-1fda-4b2a-b500-f79dc978c6c1",
        "submitter_id": "TCGA-ALCH-000017"
      },
      {
        "id": "9418f194-8741-44db-bd8f-36f4fd8c3bf2",
        "submitter_id": "TCGA-ALCH-000018"
      },
      {
        "id": "6fb2a018-c5f3-45e5-81d3-e58e7e4bf921",
        "submitter_id": "TCGA-ALCH-000019"
      },
      {
        "id": "70236972-e796-414a-9b7a-3b29b849ba7c",
        "submitter_id": "TCGA-ALCH-000020"
      },
      {
        "id": "6f78e86f-9e31-4af5-a0d9-b8970ece476d",
        "submitter_id": "TCGA-ALCH-000021"
      },
      {
        "id": "c6fcb2f0-c6bb-4b40-a761-bae3e63869cb",
        "submitter_id": "TCGA-ALCH-000002"
      },
      {
        "id": "67782964-0065-491d-b051-2ae404bb734d",
        "submitter_id": "TCGA-ALCH-000001"
      },
      {
        "id": "b45d2891-ba81-4ecc-a250-c58060934227",
        "submitter_id": "TCGA-ALCH-000012"
      }
    ]
  }
}
```

### Additional Examples

#### Example

GraphQL query for any one case in 'TCGA-LUAD' without Diagnosis information

```bare_GraphQL
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

```bare_GraphQL
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

Query for the `state` of aliquots belonging to case with `submitter_id: "TCGA-ALCH-000001"`

```bare_GraphQL
{
  aliquot(with_path_to: {type: "case", submitter_id:"TCGA-ALCH-000001"}) {
    id release_state
  }
}
```
```Response
{
  "data": {
    "aliquot": [
      {
        "id": "7af58da0-cb3e-43e2-a074-4bd8f27565ba",
        "state": "validated"
      }
    ]
  }
}
```

#### Example

GraphQL query that uses a GraphQL fragment to get specific properties from two portions and give them aliases in the response.

```bare_GraphQL
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

```bare_GraphQL
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
