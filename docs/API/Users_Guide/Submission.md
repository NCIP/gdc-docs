# Submission

## Overview

The GDC Submission API uses methods and endpoints that are distinct from those that drive the functionality of the GDC Data Portal. In particular, data and metadata that are in the process of being submitted can only be queried using [GraphQL](#querying-submitted-data-using-graphql), and not the methods described in [Search and Retrieval](Search_and_Retrieval.md).

This section describes the GDC API's submission functionality, including methods for submitting, deleting, updating, searching, and retrieving data and metadata.

## Submission endpoint

### Constructing the endpoint URL

The endpoint for submitting data to a specific project in the GDC is constructed as follows:

	https://api.gdc.cancer.gov/[API_version/]submission/Program.name/Project.code

where __[API_version/]__ is the optional API version component (see [Getting Started](Getting_Started.md)).

The values of `Program.name` and `Project.code` can be obtained from the project URL on the GDC Data Submission Portal:

	https://portal.gdc.cancer.gov/submission/Program.name/Project.code/dashboard

For more information about program name and project code see [The GDC Data Model section](/Data/Data_Model/GDC_Data_Model/#program-name-project-code-and-project-id).

#### Example

The following are URL examples for a project with `Program.name` "TCGA" and `Project.code` "ALCH":

* Submission Portal URL: `https://portal.gdc.cancer.gov/submission/TCGA/ALCH/dashboard`
* API submission endpoint (versioned): `https://api.gdc.cancer.gov/v0/submission/TCGA/ALCH`
* API submission endpoint (unversioned): `https://api.gdc.cancer.gov/submission/TCGA/ALCH`

## Submission Formats

### Metadata Formats

#### JSON and TSV

The GDC API accepts project metadata in JSON and TSV formats for the purpose of creating entities in the GDC Data Model. This includes clinical and biospecimen metadata such as disease name and stage, patient age, sample type, and certain details about the types of data collected. Upon successful data submission and project release, this metadata is indexed and becomes available for queries by data users via the GDC Data Portal and the GDC API. See [GDC Data Model](#gdc-data-model) (below) for information on accepted metadata elements and instructions for obtaining templates for metadata submission.

##### Content-Type Header

JSON is the default format for metadata submission. Submission API calls with JSON payloads should include the HTTP header `Content-Type: application/json`. Requests with TSV payloads must instead include the header `Content-Type: text/tsv`.

##### Binary Mode

Metadata files must be uploaded in raw, unencoded form. Binary mode should be used, if available, to ensure that file contents are not encoded by the upload tool before transmission. For example, when using the `curl` command-line tool, the `--data-binary` switch should be used instead of `--data`. The `--data-binary` switch is required for uploading TSV files.

### Data File Formats

The GDC API accepts a variety of data files after their metadata has been registered: BAM and FASTQ files, clinical and biospecimen supplements, slide images, and other file types. Supported data file formats are listed on the [GDC Data Dictionary](https://docs.gdc.cancer.gov/Data_Dictionary/viewer/#?view=table-entity-list&anchor=submittable_data_file).

## GDC Data Model

Submitters should review the [GDC Data Model documentation](../../Data/Data_Model/GDC_Data_Model.md) and the [GDC Data Dictionary](../../Data_Dictionary/index.md) before initiating submission.

### UUIDs

Submitters can assign UUIDs to all submittable entities other than those that correspond to files (entities in categories `data_file` or `metadata_file`). If the submitter does not provide a UUID, it will be assigned by the GDC and returned in the API response upon successful completion of the submission transaction. See [Appendix C](Appendix_C_Format_of_Submission_Requests_and_Responses.md) for details of the API response format. To learn more about UUIDs see the [GDC Data Model documentation](../../Data/Data_Model/GDC_Data_Model.md#uuids).

### Submitter IDs

In addition to `id`, many entities also include a `submitter_id` field. This field can contain any string (e.g. a "barcode") that the submitter wishes to use to identify the entity. Typically this string identifies a corresponding entry in submitter's records. The GDC's only requirement with respect to `submitter_id` is that it be a string that is unique for all entities within a project. The GDC Submission API requires a `submitter_id` for most entities.

>**Note:** For `case` entities, `submitter_id` must correspond to a `submitted_subject_id` of a study participant registered with the project in dbGaP.

### GDC Data Dictionary Endpoints

Information in the [GDC Data Dictionary](../../Data_Dictionary/index.md) can be accessed programmatically as described below.

#### Submission Templates

Submission templates are accessible programmatically at the `templates` endpoint. Template format (`json`, `tsv` or `csv`) is specified using the `format` parameter.

For example, the JSON template for `case` entities can be obtained from:

	https://api.gdc.cancer.gov/v0/submission/template/case?format=json

In addition to `case`, templates for the following entities can be downloaded


__Biospecimen:__
```
sample
portion
analyte
aliquot
read_group
```
__Clinical:__
```
slide
demographic
diagnosis
exposure
family_history
treatment
follow_up
molecular_test
```

__Data Files:__
```
analysis_metadata
biospecimen_supplement
clinical_supplement
experiment_metadata
pathology_report
run_metadata
slide_image
submitted_unaligned_reads
submitted_aligned_reads
submitted_genomic_profile
```

#### Entity JSON Schemas

The entire collection of GDC entity schemas can be downloaded from the `dictionary` endpoint:

	https://api.gdc.cancer.gov/v0/submission/_dictionary/_all

Individual schemas can be downloaded by specifying entity type. For example, the JSON schema for `case` entities can be found at:

	https://api.gdc.cancer.gov/v0/submission/_dictionary/case


## Making Requests to the Submission API

Requests to create or update entities in the GDC must specify the entity `type`, the entity `id` or `submitter_id`, relationships (links) that the entity has to existing entities in the [GDC Data Model](../../Data/Data_Model/GDC_Data_Model.md), and entity properties as defined by the [GDC Data Dictionary](../../Data_Dictionary/index.md). To delete entities, only the `id` property is required. The general format of GDC API submission requests and responses is provided in [Appendix C](Appendix_C_Format_of_Submission_Requests_and_Responses.md).

## Submission Transactions

Submission of data to the GDC involves a series of transactions initiated by the submitter, that create and link entities according to the [GDC Data Model](../../Data/Data_Model/GDC_Data_Model.md). With the exception of `program`, which is an administrative entity created by the GDC, all new entities must be linked, at creation, to existing entities or to new entities being created in the same transaction. For example, a submitter cannot create a `portion` entity unless the submitter either (1) has previously created the corresponding `case` and `sample` entities, or (2) is creating those entities in the same transaction. This also means that entities cannot be deleted if they have "child" entities attached to them.

If multiple entities are being created and/or updated in a transaction, and an error is encountered for one of the entities, then the transaction will fail and no changes will be made to the GDC.

### Dry Run Transactions

The `submission` endpoint provides a `_dry_run` mode that simulates submission transactions without making changes to the GDC. This mode is activated by appending `/_dry_run` to the end of a submission endpoint.

The following is an example of a POST request, that simulates creating an entity in dry run mode:

=== "Request"

    ```json
    {
      "project_id": "GDC-INTERNAL",
      "type": "case",
      "submitter_id": "GDC-INTERNAL-000093",
      "disease_type": "Blood Vessel Tumors",
      "primary_site": "Base of tongue",
       "projects": {
         "code": "INTERNAL"
      }
    }

    ```

=== "Command"

    ```shell
    token=$(<gdc-token-text-file.txt)

    curl --header "X-Auth-Token: $token" --request POST --data-binary @Request --header 'Content-Type: application/json' https://api.gdc.cancer.gov/v0/submission/GDC/INTERNAL/_dry_run
    ```

=== "Response"

    ```json
    {
      "cases_related_to_created_entities_count": 0,
      "cases_related_to_updated_entities_count": 0,
      "code": 200,
      "created_entity_count": 1,
      "entities": [
        {
          "action": "create",
          "errors": [],
          "id": "bfc1fb29-28db-4137-8379-3d1693ce3423",
          "related_cases": [],
          "type": "case",
          "unique_keys": [
            {
              "project_id": "GDC-INTERNAL",
              "submitter_id": "GDC-INTERNAL-000093"
            }
          ],
          "valid": true,
          "warnings": []
        }
      ],
      "entity_error_count": 0,
      "message": "Transaction would have been successful. User selected dry run option, transaction aborted, no data written to database.",
      "success": true,
      "transaction_id": 5834800,
      "transactional_error_count": 0,
      "transactional_errors": [],
      "updated_entity_count": 0
    }
    ```

#### Dry Run Commit

For convenience, the GDC enables users to commit earlier `_dry_run` transactions instead of uploading the same data again to execute the changes. This `commit` action is allowed on transactions that (1) have not been previously committed and (2) were successful `dry_run` transactions.

Note that the `commit` action is a separate transaction with its own transaction id, and it can be executed [asynchronously](#asynchronous-transactions). If the state of the submission project has changed in a way that would make the original `_dry_run` transaction invalid if it were run again (e.g. entities with the same `submitter_id` have since been created in another transaction), then the `commit` action will fail.

To commit a transaction, submit a POST or PUT request to `/submission/Program.name/Project.code/transactions/transaction_id/commit`, replacing `Program.name`, `Project.code`, and `transaction_id` with values associated with the transaction.

=== "Command"

    ```shell
    token=$(<gdc-token-text-file.txt)

    curl --header "X-Auth-Token: $token" --request POST https://api.gdc.cancer.gov/v0/submission/GDC/INTERNAL/transactions/467/commit?async=true
    ```

=== "Response"

    ```json
    {
      "code": 200,
      "message": "Transaction submitted.",
      "transaction_id": 468,
    }
    ```



#### Dry Run Close

The GDC Submission API also provides a `close` action on `_dry_run` transactions. This `close` action is allowed on `_dry_run` transactions that have not been previously closed. Closing a `_dry_run` transaction prevents it from being committed in the future.

To close a transaction, submit a POST or PUT request to `/submission/Program.name/Project.code/transactions/transaction_id/close`, replacing `Program.name`, `Project.code`, and `transaction_id` with values associated with the transaction.

=== "Command"

    ```shell
    token=$(<gdc-token-text-file.txt)

    curl --header "X-Auth-Token: $token" --request POST https://api.gdc.cancer.gov/v0/submission/GDC/INTERNAL/transactions/467/close
    ```

=== "Response"

    ```json
    {
        "code": 200,
        "message": "Closed transaction.",
        "transaction_id": 467
    }
    ```


### Asynchronous Transactions

The `submission` endpoint provides an asynchronous mode that provides immediate response and executes submission transactions in the background. This mode is activated by appending `?async=true` to the end of a submission endpoint.  The API will respond with the `transaction_id` which can be used to look up the result of the transaction at a later time via the [GraphQL](#querying-submitted-data-using-graphql) endpoint.  If the server has too many asynchronous jobs scheduled already, your request to schedule a transaction may fail.

#### Example

The following is an example of a PUT request, that creates a case asynchronously:

=== "Request"

    ```json
    {
      "project_id": "GDC-INTERNAL",
      "type": "case",
      "submitter_id": "GDC-INTERNAL-000093",
      "disease_type": "Blood Vessel Tumors",
      "primary_site": "Base of tongue",
       "projects": {
         "code": "INTERNAL"
      }
    }
    ```

=== "Command"

    ```shell
    token=$(<gdc-token-text-file.txt)

    curl --header "X-Auth-Token: $token" --request POST --data-binary @Request --header 'Content-Type: application/json' https://api.gdc.cancer.gov/v0/submission/GDC/INTERNAL?async=true
    ```

=== "Response"

    ```json
    {
      "code": 200,
      "message": "Transaction submitted.",
      "transaction_id": 467,
    }
    ```

The following is a [GraphQL](#querying-submitted-data-using-graphql) request that looks up the state of the above transaction:

=== "GraphQL_Request"

    ```GraphQL
    query {
      transaction_log(id: 467) {
        is_dry_run
        committed_by
        state
      }
    }
    ```

=== "GraphQL_Response"

    ```GraphQL
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

The following transaction fields can be queried using [GraphQL](#querying-submitted-data-using-graphql) and are helpful in determining the status of a transaction:

|Field|Type|Description|
|---|---|---|
|`id`|ID|Transaction identifier|
|`is_dry_run`|Boolean|Indicates whether the transaction is a dry run|
|`closed`|Boolean|For dry run transactions, indicates whether the transaction has been closed to prevent it from being committed in the future.|
|`committable`|Boolean|Indicates whether the transaction can be committed (i.e. it is a successful dry run transaction that has not been committed previously and has not been closed)|
|`state`|String|Indicates the state of the transaction: `PENDING`, `SUCCEEDED`, `FAILED` (due to user error), or `ERRORED` (due to system error)|
|`committed_by`|ID|The ID of the transaction that committed this transaction|

>**Note:** To check whether a dry run transaction was committed successfully, check the `state` of the transaction that executed the commit. The `state` of the dry run transaction itself does not represent the status of a subsequent commit.

## Creating and Updating Entities

The GDC Submission API supports HTTP POST and HTTP PUT methods for creating entities:

* **POST** will create entities that do not exist, and will fail if any of the entities in the transaction already exist in the GDC.

* **PUT** will create new entities and update existing entities, and identify which entities were created or updated in the API response.

The GDC suggests using POST for creating new entities, and using PUT only for updating entities. This helps to avoid inadvertent entity updates that can occur when using PUT for creating entities.

>**Note:** Once a relationship has been created between two entities, it cannot be removed by updating an entity. To remove a relationship, the child entity must be [deleted](#deleting-entities).


### Example: Creating and Updating Case Entities (JSON)

In this example, a case entity is created using POST. Then an attempt is made to create the same entity again using POST, resulting in an error. Then the originally created entity is updated (with the same information) using PUT.

The JSON in the request was generated using the `case` JSON template that can be obtained from the [GDC Data Dictionary Viewer](../../Data_Dictionary/index.md) and from `https://api.gdc.cancer.gov/v0/submission/template/case?format=json`.

>**Note:** For `case` entities, `submitter_id` must correspond to a `submitted_subject_id` of a study participant registered with the project in dbGaP.


=== "Request1"

    ```json
    {
      "project_id": "GDC-INTERNAL",
      "type": "case",
      "submitter_id": "GDC-INTERNAL-000093",
      "disease_type": "Blood Vessel Tumors",
      "primary_site": "Base of tongue",
       "projects": {
         "code": "INTERNAL"
      }
    }
    ```

=== "Command1"

    ```shell
    token=$(<gdc-token-text-file.txt)

    curl --header "X-Auth-Token: $token" --request POST --data-binary @Request --header 'Content-Type: application/json' https://api.gdc.cancer.gov/v0/submission/GDC/INTERNAL
    ```

=== "Response1"

    ```json
    {
      "cases_related_to_created_entities_count": 0,
      "cases_related_to_updated_entities_count": 0,
      "code": 200,
      "created_entity_count": 1,
      "entities": [
        {
          "action": "create",
          "errors": [],
          "id": "bfc1fb29-28db-4137-8379-3d1693ce3423",
          "related_cases": [],
          "type": "case",
          "unique_keys": [
            {
              "project_id": "GDC-INTERNAL",
              "submitter_id": "GDC-INTERNAL-000093"
            }
          ],
          "valid": true,
          "warnings": []
        }
      ],
      "entity_error_count": 0,
      "message": "Transaction would have been successful. User selected dry run option, transaction aborted, no data written to database.",
      "success": true,
      "transaction_id": 5834800,
      "transactional_error_count": 0,
      "transactional_errors": [],
      "updated_entity_count": 0
    }
    ```

=== "Command2"

    ```shell
    curl --header "X-Auth-Token: $token" --request POST --data-binary @Request --header 'Content-Type: application/json' https://api.gdc.cancer.gov/v0/submission/TCGA/ALCH
    ```

=== "Response2"

    ```json
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
              "message": "Cannot create an entity with an id that already exists.",
              "type": "NOT_UNIQUE"
            }
          ],
          "id": null,
          "related_cases": [],
          "type": "case",
          "unique_keys": [
            {
              "project_id": "GDC-INTERNAL",
              "submitter_id": "GDC-INTERNAL-000093"
            }
          ],
          "valid": false,
          "warnings": []
        }
      ],
      "entity_error_count": 1,
      "message": "Transaction aborted due to 1 invalid entity.",
      "success": false,
      "transaction_id": 5834802,
      "transactional_error_count": 0,
      "transactional_errors": [],
      "updated_entity_count": 0
    }
    ```

=== "Command3"

    ```shell
    curl --header "X-Auth-Token: $token" --request PUT --data-binary @Request --header 'Content-Type: application/json' https://api.gdc.cancer.gov/v0/submission/GDC/INTERNAL
    ```

=== "Response3"

    ```json
    {
      "cases_related_to_created_entities_count": 0,
      "cases_related_to_updated_entities_count": 0,
      "code": 200,
      "created_entity_count": 0,
      "entities": [
        {
          "action": "update",
          "errors": [],
          "id": "a35b8e26-3b43-4203-9d33-44c2b351f177",
          "related_cases": [],
          "type": "case",
          "unique_keys": [
            {
              "project_id": "GDC-INTERNAL",
              "submitter_id": "GDC-INTERNAL-000093"
            }
          ],
          "valid": true,
          "warnings": []
        }
      ],
      "entity_error_count": 0,
      "message": "Transaction successful.",
      "success": true,
      "transaction_id": 5834803,
      "transactional_error_count": 0,
      "transactional_errors": [],
      "updated_entity_count": 1
    }
    ```



### Example: Creating an Aliquot Entity (JSON)

In this example, an `aliquot` entity and a `sample` entity are created in a single transaction. The `aliquot` is linked to `sample` which is linked to `case`. The first request is an example of using `submitter_id` properties to link entities together. The second request is an example of using UUIDs for creating the links.

#### Request 1: Creating Links Using submitter_id

=== "Request"

    ```json
    [
      {
        "type": "sample",
        "submitter_id": "GDC-INTERNAL-000093-SAMPLE000093",
        "tissue_type": "Tumor",
        "preservation_method": "Fresh",
        "specimen_type": "Whole Bone Marrow",
        "tumor_descriptor": "Primary",
        "cases": {
          "submitter_id": "GDC-INTERNAL-000093"
        }
      },
      {
        "type": "aliquot",
        "submitter_id": "GDC-INTERNAL-000093-SAMPLE000093-ALIQUOT000093",
        "samples": {
          "submitter_id": "GDC-INTERNAL-000093-SAMPLE000093"
        }
      }
    ]
    ```

=== "Command"

    ```shell
    token=$(<gdc-token-text-file.txt)

    curl --header "X-Auth-Token: $token" --request POST --data-binary @Request --header 'Content-Type: application/json' https://api.gdc.cancer.gov/v0/submission/GDC/INTERNAL
    ```

=== "Response"

    ```json
    {
      "cases_related_to_created_entities_count": 1,
      "cases_related_to_updated_entities_count": 0,
      "code": 201,
      "created_entity_count": 2,
      "entities": [
        {
          "action": "create",
          "errors": [],
          "id": "0a877533-0c85-4a7e-9309-733ccf295c1b",
          "related_cases": [
            {
              "id": "a00f076e-d694-47dd-8e50-24c28e90fd6a",
              "submitter_id": "GDC-INTERNAL-000093"
            }
          ],
          "type": "sample",
          "unique_keys": [
            {
              "project_id": "GDC-INTERNAL",
              "submitter_id": "GDC-INTERNAL-000093-SAMPLE000093"
            }
          ],
          "valid": true,
          "warnings": []
        },
        {
          "action": "create",
          "errors": [],
          "id": "45c57067-c92d-453b-8b6d-14a3fe08f802",
          "related_cases": [
            {
              "id": "a00f076e-d694-47dd-8e50-24c28e90fd6a",
              "submitter_id": "GDC-INTERNAL-000093"
            }
          ],
          "type": "aliquot",
          "unique_keys": [
            {
              "project_id": "GDC-INTERNAL",
              "submitter_id": "GDC-INTERNAL-000093-SAMPLE000093-ALIQUOT000093"
            }
          ],
          "valid": true,
          "warnings": []
        }
      ],
      "entity_error_count": 0,
      "message": "Transaction successful.",
      "success": true,
      "transaction_id": 5835160,
      "transactional_error_count": 0,
      "transactional_errors": [],
      "updated_entity_count": 0
    }
    ```

#### Request 2: Creating Links Using UUID

=== "Request"

    ```json
    [
      {
        "type": "sample",
        "submitter_id": "GDC-INTERNAL-000093-SAMPLE000093",
        "tissue_type": "Tumor",
        "preservation_method": "Fresh",
        "specimen_type": "Whole Bone Marrow",
        "tumor_descriptor": "Primary",
        "cases": {
          "id": "a00f076e-d694-47dd-8e50-24c28e90fd6a"
        }
      },
      {
        "type": "aliquot",
        "submitter_id": "GDC-INTERNAL-000093-SAMPLE000093-ALIQUOT000093",
        "samples": {
          "submitter_id": "GDC-INTERNAL-000093-SAMPLE000093"
        }
      }
    ]
    ```

=== "Command"

    ```shell
    token=$(<gdc-token-text-file.txt)

    curl --header "X-Auth-Token: $token" --request POST --data-binary @Request --header 'Content-Type: application/json' https://api.gdc.cancer.gov/submission/GDC/INTERNAL
    ```

=== "Response"

    ```json
    {
      "cases_related_to_created_entities_count": 1,
      "cases_related_to_updated_entities_count": 0,
      "code": 201,
      "created_entity_count": 2,
      "entities": [
        {
          "action": "create",
          "errors": [],
          "id": "9684fd7c-97b5-42a2-b350-2d86d41bbfdb",
          "related_cases": [
            {
              "id": "a00f076e-d694-47dd-8e50-24c28e90fd6a",
              "submitter_id": "GDC-INTERNAL-000093"
            }
          ],
          "type": "sample",
          "unique_keys": [
            {
              "project_id": "GDC-INTERNAL",
              "submitter_id": "GDC-INTERNAL-000093-SAMPLE000093"
            }
          ],
          "valid": true,
          "warnings": []
        },
        {
          "action": "create",
          "errors": [],
          "id": "cd5613ef-acb8-4b56-af4b-8bf3ab0e09d8",
          "related_cases": [
            {
              "id": "a00f076e-d694-47dd-8e50-24c28e90fd6a",
              "submitter_id": "GDC-INTERNAL-000093"
            }
          ],
          "type": "aliquot",
          "unique_keys": [
            {
              "project_id": "GDC-INTERNAL",
              "submitter_id": "GDC-INTERNAL-000093-SAMPLE000093-ALIQUOT000093"
            }
          ],
          "valid": true,
          "warnings": []
        }
      ],
      "entity_error_count": 0,
      "message": "Transaction successful.",
      "success": true,
      "transaction_id": 5835208,
      "transactional_error_count": 0,
      "transactional_errors": [],
      "updated_entity_count": 0
    }
    ```

### Example: Creating Two Samples (TSV)

In this example, a TSV file containing metadata for two samples is uploaded to the GDC in dry run mode.

=== "Request"

    ```
    type	project_id	submitter_id	cases.submitter_id	specimen_type	tissue_type	tumor_descriptor	preservation_method
    sample	GDC-INTERNAL	GDC-INTERNAL-000093-sampleA	GDC-INTERNAL-000093	Solid Tissue	Tumor	Primary	Frozen
    sample	GDC-INTERNAL	GDC-INTERNAL-000093-sampleB	GDC-INTERNAL-000093	Solid Tissue	Normal	Not Reported	Frozen
    ```

=== "Command"

    ```shell
    curl --header "X-Auth-Token: $token" --header 'Content-Type: text/tsv' --request PUT --data-binary @Samples.tsv 'https://api.gdc.cancer.gov/submission/GDC/INTERNAL/_dry_run'
    ```

=== "Reponse"

    ```json
    {
      "cases_related_to_created_entities_count": 1,
      "cases_related_to_updated_entities_count": 0,
      "code": 200,
      "created_entity_count": 2,
      "entities": [
        {
          "action": "create",
          "errors": [],
          "id": "c3d401f1-d505-4240-b801-dc2b389ddea1",
          "related_cases": [
            {
              "id": "a00f076e-d694-47dd-8e50-24c28e90fd6a",
              "submitter_id": "GDC-INTERNAL-000093"
            }
          ],
          "type": "sample",
          "unique_keys": [
            {
              "project_id": "GDC-INTERNAL",
              "submitter_id": "GDC-INTERNAL-000093-sampleA"
            }
          ],
          "valid": true,
          "warnings": []
        },
        {
          "action": "create",
          "errors": [],
          "id": "d3c4be95-9c69-4e8e-9f37-61db455ded7a",
          "related_cases": [
            {
              "id": "a00f076e-d694-47dd-8e50-24c28e90fd6a",
              "submitter_id": "GDC-INTERNAL-000093"
            }
          ],
          "type": "sample",
          "unique_keys": [
            {
              "project_id": "GDC-INTERNAL",
              "submitter_id": "GDC-INTERNAL-000093-sampleB"
            }
          ],
          "valid": true,
          "warnings": []
        }
      ],
      "entity_error_count": 0,
      "message": "Transaction would have been successful. User selected dry run option, transaction aborted, no data written to database.",
      "success": true,
      "transaction_id": 5835321,
      "transactional_error_count": 0,
      "transactional_errors": [],
      "updated_entity_count": 0
    }
    ```

### Example: Updating a Sample Entity (JSON)

Entities can be updated using a very similar process to what is shown above.  

#### Updating a sample

New nodes are created in Request1.  Nodes in state `validated` are updated in Request2.

=== "Request1"

    ```json
    [
       {
        "type": "case",
        "submitter_id": "QA-REGRESSION-0002",  
        "projects": {
        "code": "REGRESSION"
      }
        },
        {
        "type": "sample",
        "submitter_id": "QA-REGRESSION-0002-SAMPLE000001",
        "sample_type": "Primary Tumor",
        "sample_type_id": "01",
        "cases": {
          "submitter_id": "QA-REGRESSION-0002"
        }
      },
      {
        "type": "aliquot",
        "submitter_id": "QA-REGRESSION-0002-SAMPLE000001-ALIQUOT000001",
        "samples": {
          "submitter_id": "QA-REGRESSION-0002-SAMPLE000001"
        }
      }
    ]
    ```

=== "Conmmand1"

    ```shell
    token=$(<gdc-token-text-file.txt)

    curl --header "X-Auth-Token: $token" --request POST --data-binary @sample.json --header 'Content-Type: application/json' https://api.gdc.cancer.gov/v0/submission/QA/REGRESSION
    ```

=== "Response1"

    ```json
    {
      "cases_related_to_created_entities_count": 1,
      "cases_related_to_updated_entities_count": 0,
      "code": 201,
      "created_entity_count": 3,
      "entities": [
        {
          "action": "create",
          "errors": [],
          "id": "3a750ae8-8e63-472e-852e-8e514a0c1550",
          "related_cases": [],
          "type": "case",
          "unique_keys": [
            {
              "project_id": "QA-REGRESSION",
              "submitter_id": "QA-REGRESSION-0002"
            }
          ],
          "valid": true,
          "warnings": []
        },
        {
          "action": "create",
          "errors": [],
          "id": "8a1872e6-c5e6-4f39-b9fe-15ecf45715c7",
          "related_cases": [
            {
              "id": "3a750ae8-8e63-472e-852e-8e514a0c1550",
              "submitter_id": "QA-REGRESSION-0002"
            }
          ],
          "type": "sample",
          "unique_keys": [
            {
              "project_id": "QA-REGRESSION",
              "submitter_id": "QA-REGRESSION-0002-SAMPLE000001"
            }
          ],
          "valid": true,
          "warnings": []
        },
        {
          "action": "create",
          "errors": [],
          "id": "e9279137-92b4-41ab-be28-a03e32e6fac7",
          "related_cases": [
            {
              "id": "3a750ae8-8e63-472e-852e-8e514a0c1550",
              "submitter_id": "QA-REGRESSION-0002"
            }
          ],
          "type": "aliquot",
          "unique_keys": [
            {
              "project_id": "QA-REGRESSION",
              "submitter_id": "QA-REGRESSION-0002-SAMPLE000001-ALIQUOT000001"
            }
          ],
          "valid": true,
          "warnings": []
        }
      ],
      "entity_error_count": 0,
      "message": "Transaction successful.",
      "success": true,
      "transaction_id": 920117,
      "transactional_error_count": 0,
      "transactional_errors": [],
      "updated_entity_count": 0
    }
    ```

=== "Request2"

    ```json
    [
       {
        "type": "case",
        "submitter_id": "QA-REGRESSION-0002",  
        "projects": {
        "code": "REGRESSION"
      }
        },
        {
        "type": "sample",
        "submitter_id": "QA-REGRESSION-0002-SAMPLE000001",
        "sample_type": "Primary Tumor",
        "days_to_collection":5,
        "sample_type_id": "01",
        "cases": {
          "submitter_id": "QA-REGRESSION-0002"
        }
      },
      {
        "type": "aliquot",diff
        "submitter_id": "QA-REGRESSION-0002-SAMPLE000001-ALIQUOT000001",
        "samples": {
          "submitter_id": "QA-REGRESSION-0002-SAMPLE000001"
        }
      }
    ]
    ```

=== "Command2"

    ```shell
    token=$(<gdc-token-text-file.txt)

    curl --header "X-Auth-Token: $token" --request PUT --data-binary @sample2.json --header 'Content-Type: application/json' https://api.gdc.cancer.gov/v0/submission/QA/REGRESSION
    ```

=== "Response2"

    ```json
    {
      "cases_related_to_created_entities_count": 0,
      "cases_related_to_updated_entities_count": 1,
      "code": 200,
      "created_entity_count": 0,
      "entities": [
        {
          "action": "update",
          "errors": [],
          "id": "3a750ae8-8e63-472e-852e-8e514a0c1550",
          "related_cases": [],
          "type": "case",
          "unique_keys": [
            {
              "project_id": "QA-REGRESSION",
              "submitter_id": "QA-REGRESSION-0002"
            }
          ],
          "valid": true,
          "warnings": []
        },
        {
          "action": "update",
          "errors": [],
          "id": "8a1872e6-c5e6-4f39-b9fe-15ecf45715c7",
          "related_cases": [
            {
              "id": "3a750ae8-8e63-472e-852e-8e514a0c1550",
              "submitter_id": "QA-REGRESSION-0002"
            }
          ],
          "type": "sample",
          "unique_keys": [
            {
              "project_id": "QA-REGRESSION",
              "submitter_id": "QA-REGRESSION-0002-SAMPLE000001"
            }
          ],
          "valid": true,
          "warnings": []
        },
        {
          "action": "update",
          "errors": [],
          "id": "e9279137-92b4-41ab-be28-a03e32e6fac7",
          "related_cases": [
            {
              "id": "3a750ae8-8e63-472e-852e-8e514a0c1550",
              "submitter_id": "QA-REGRESSION-0002"
            }
          ],
          "type": "aliquot",
          "unique_keys": [
            {
              "project_id": "QA-REGRESSION",
              "submitter_id": "QA-REGRESSION-0002-SAMPLE000001-ALIQUOT000001"
            }
          ],
          "valid": true,
          "warnings": []
        }
      ],
      "entity_error_count": 0,
      "message": "Transaction successful.",
      "success": true,
      "transaction_id": 920120,
      "transactional_error_count": 0,
      "transactional_errors": [],
      "updated_entity_count": 3
    }
    ```

## Retrieving Entities

### Entities Endpoint

JSON objects representing submitted entities can be retrieved with a GET request to the `entities` endpoint. This endpoint retrieves entities by UUID. A single UUID or a comma-separated list of UUIDs can be passed to this endpoint as a query.

=== "Command"

    ```shell
    token=$(<gdc-token-text-file.txt)

    curl --header "X-Auth-Token: $token" https://api.gdc.cancer.gov/v0/submission/TCGA/ALCH/entities/fbf69646-5904-4f95-92d6-692bde658f05
    ```

=== "Response"

    ```json
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

The `export` endpoint provides additional functionality for exporting entities from the GDC submission system. The `ids` parameter accepts a UUID or a comma-separated list of UUIDs. The `format` parameter allows the user to specify the preferred format of the API response: JSON, TSV, or CSV. When the `with_children` parameter is set to `with_children`, the response includes the metadata stored in all "child" entities of the entity being requested. The `export` endpoint accepts GET requests.

=== "Command"

    ```shell
    token=$(<gdc-token-text-file.txt)


    curl --header "X-Auth-Token: $token" 'https://api.gdc.cancer.gov/v0/submission/TCGA/ALCH/export?ids=11f8321-832f-4a8b-8384-a2f6256557e0&format=json&with_children=with_children'
    ```

=== "Response"

    ```json
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

## Patching Entitites

The GDC Submission API supports the HTTP PATCH method for updating existing entities with additional fields.

**PATCH** can be used to add extra fields to an existing entity, without requiring the submission of required fields.

The PATCH method cannot be used to create new entities, and the provided submitter_id must match an existing submitter_id.

#### Example: Creating a new demographic entity using POST

=== "Request1"

    ```json
    {
      "type": "demographic",
      "submitter_id": "demographic7892",
      "cases": {
        "submitter_id": "GDC-INTERNAL-000073"
      },
      "ethnicity": "not reported",
      "sex_at_birth": "male",
      "race": "white",
      "vital_status": "Dead"
    }
    ```

=== "Command1"

    ```shell
    token=$(<gdc-token-text-file.txt)

    curl --header "X-Auth-Token: $token" --request POST --data-binary @Request --header 'Content-Type: application/json' https://api.gdc.cancer.gov/v0/submission/GDC/INTERNAL
    ```

=== "Response1"

    ```json
    {
      "cases_related_to_created_entities_count": 1,
      "cases_related_to_updated_entities_count": 0,
      "code": 201,
      "created_entity_count": 1,
      "entities": [
        {
          "action": "create",
          "errors": [],
          "id": "4e4f29a3-5325-47ef-a583-a251677ed29a",
          "related_cases": [
            {
              "id":"71d17c1f-8985-4b2f-bb63-1c39cb6562d5",
              "submitter_id":"GDC-INTERNAL-000073"
            }
          ],
          "type": "demographic",
          "unique_keys": [
            {
              "project_id": "GDC-INTERNAL",
              "submitter_id": "demographic7892"
            }
          ],
          "valid": true,
          "warnings": []
        }
      ],
      "entity_error_count": 0,
      "message": "Transaction successful.",
      "success": true,
      "transaction_id": 6357750,
      "transactional_error_count": 0,
      "transactional_errors": [],
      "updated_entity_count": 0
    }
    ```

#### Example: Updating the existing demographic entity using PATCH

=== "Request2"

    ```json
    {
      "type": "demographic",
      "submitter_id": "demographic7892",
      "cause_of_death": "Infection",
      "cause_of_death_source": "Death Certificate",
      "country_of_birth": "Antigua and Barbuda",
      "country_of_residence_at_enrollment": "Antigua and Barbuda"
    }
    ```

=== "Command2"

    ```shell
    token=$(<gdc-token-text-file.txt)

    curl --header "X-Auth-Token: $token" --request PATCH --data-binary @Request --header 'Content-Type: application/json' https://api.gdc.cancer.gov/v0/submission/GDC/INTERNAL
    ```

=== "Response2"

    ```json
    {
      "cases_related_to_created_entities_count": 0,
      "cases_related_to_updated_entities_count": 1,
      "code": 200,
      "created_entity_count": 0,
      "entities": [
        {
          "action": "update",
          "errors": [],
          "id": "4e4f29a3-5325-47ef-a583-a251677ed29a",
          "related_cases": [
            {
              "id":"71d17c1f-8985-4b2f-bb63-1c39cb6562d5",
              "submitter_id":"GDC-INTERNAL-000073"
            }
          ],
          "type": "demographic",
          "unique_keys": [
            {
              "project_id": "GDC-INTERNAL",
              "submitter_id": "demographic7892"
            }
          ],
          "valid": true,
          "warnings": []
        }
      ],
      "entity_error_count": 0,
      "message": "Transaction successful.",
      "success": true,
      "transaction_id": 6357751,
      "transactional_error_count": 0,
      "transactional_errors": [],
      "updated_entity_count": 1
    }
    ```


## Deleting Entities

The `entities` endpoint can also be used to delete entities. This is accomplished using a DELETE request to the endpoint, specifying the entity's UUID. If an entity cannot be deleted because it is linked to child entities, the GDC Submission API will respond with an error providing a list of entities that must be deleted prior to deleting the subject entity.

A subgraph (a parent along with all of its child entities) can be deleted in a single transaction by passing a comma-separated list of UUIDs to the `entities` endpoint.

Entities in submitted state (assigned when the project has been submitted) cannot be deleted.

=== "Shell"

    ```Shell
    token=$(<gdc-token-text-file.txt)

    curl --header "X-Auth-Token: $token" --request DELETE https://api.gdc.cancer.gov/v0/submission/TCGA/ALCH/entities/67782964-0065-491d-b051-2ae404bb734d
    ```

=== "Response"

    ```json
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

Experimental data files like BAM and FASTQ can be uploaded directly to the API using the `files` endpoint, by specifying the UUID of the corresponding `data_file` entity. Binary upload mode must be used if available. Uploading large files may be more efficiently performed using the [GDC Data Transfer Tool](/Data_Transfer_Tool/Users_Guide/Getting_Started/).

```
token=$(<gdc-token-text-file.txt)

curl --header "X-Auth-Token: $token" --output needed_to_show_progress_bar.log --request PUT --data-binary @GDC-INTERNAL-000084-S1-Q1-RG1.fastq.zip https://api.gdc.cancer.gov/v0/submission/GDC/INTERNAL/files/c414a205-376e-4993-af48-2a4689eb433e && rm needed_to_show_progress_bar.log

	# "&& rm needed_to_show_progress_bar.log" at the end of the command above
	# removes the temporary file required to show upload progress bar. This
	# will not work on Windows platforms. Windows users must remove this
	# string and can delete the file manually.
```

#### Upload Manifest

The `manifest` endpoint generates a manifest for uploading files using the GDC Data Transfer Tool. It requires a comma-separated list of file UUIDs to generate a manifest.

```
https://api.gdc.cancer.gov/v0/submission/PROGRAM/PROJECT/manifest?ids=bf0751ca-fc3b-4760-b876-0fefce040be5,90163202-cfd7-4f6a-8214-e7e4e924d3a6
```
### Uploading New Versions of Data Files

The GDC Submission system supports submitting updated versions of files.  For example, you may want to submit an updated version of a clinical supplement file that contains new clinical information about a patient.  If a file is in file_state `validated` then you would simply delete and upload a new copy of this file.  No additional version of the file will be created in this case.  The UUID of the node stays the same.

However, if a file is in file_state `submitted` or `validated` and state `released` then a different process is required.  In this situation simply upload a new template containing updated metadata (e.g. md5sum or file_size).  A new node (with a new UUID) will automatically be created that is linked to the previous version.  Once this new file is indexed and released to users they will be able to query the new UUID in the /files endpoint and both versions' UUID in the files/versions or /history endpoint.  In the example below we register a file, upload the file, and register a new version of this file.

=== "Request"

    ```json
    [
      {
        "data_type": "Clinical Supplement",
        "file_name": "nationwidechildrens.org_clinical.TCGA-4G-AAZT-.xml",
        "md5sum": "ecaaa87613ba03c971bfefdb6f693959",
        "data_format": "BCR XML",
        "submitter_id": "nationwidechildrens.org_CHOL.bio.Level_1.428.25.0.tar.gz_nationwidechildrens.org_clinical.TCGA-4G-AAZT-.xml",
        "archives": [],
        "data_category": "Clinical",
        "file_size": 39195,
        "cases": [
          {
            "id": "b10c64c2-7fd2-4210-b975-034affb14b57",
            "submitter_id": "TCGA-4G-AAZT"
          }
        ],
        "project_id": "TCGA-CHOL",
        "type": "clinical_supplement"
      }
    ]
    ```

=== "Command"

    ```shell
    curl --header "X-Auth-Token: $token" --header 'Content-Type: json' --request PUT --data-binary @clin.json 'https://api.gdc.cancer.gov/submission/TCGA/CHOL'
    ```

=== "Response"

    ```json
    {
      "cases_related_to_created_entities_count": 1,
      "cases_related_to_updated_entities_count": 0,
      "code": 200,
      "created_entity_count": 1,
      "entities": [
        {
          "action": "create",
          "errors": [],
          "id": "d65c15d9-9e33-4a0b-863d-605ad6155506",
          "related_cases": [
            {
              "id": "b10c64c2-7fd2-4210-b975-034affb14b57",
              "submitter_id": "TCGA-4G-AAZT"
            }
          ],
          "type": "clinical_supplement",
          "unique_keys": [
            {
              "project_id": "TCGA-CHOL",
              "submitter_id": "nationwidechildrens.org_CHOL.bio.Level_1.428.25.0.tar.gz_nationwidechildrens.org_clinical.TCGA-4G-AAZT-.xml"
            }
          ],
          "valid": true,
          "warnings": []
        }
      ],
      "entity_error_count": 0,
      "message": "Transaction successful.",
      "success": true,
      "transaction_id": 922606,
      "transactional_error_count": 0,
      "transactional_errors": [],
      "updated_entity_count": 0
    }
    ```

=== "Request2"

    ```json
    [
      {
        "data_type": "Clinical Supplement",
        "file_name": "nationwidechildrens.org_clinical.TCGA-4G-AAZT-.xml",
        "md5sum": "93e306e5e621d3cacb363e2be96ca3cd",
        "data_format": "BCR XML",
        "submitter_id": "nationwidechildrens.org_CHOL.bio.Level_1.428.25.0.tar.gz_nationwidechildrens.org_clinical.TCGA-4G-AAZT-.xml",
        "archives": [],
        "data_category": "Clinical",
        "file_size": 39197,
        "cases": [
          {
            "id": "b10c64c2-7fd2-4210-b975-034affb14b57",
            "submitter_id": "TCGA-4G-AAZT"
          }
        ],
        "project_id": "TCGA-CHOL",
        "type": "clinical_supplement"
      }
    ]
    ```

=== "Command2"

    ```shell
    curl --header "X-Auth-Token: $token" --header 'Content-Type: json' --request PUT --data-binary @clin_v2.json 'https://api.gdc.cancer.gov/submission/TCGA/CHOL'
    ```

=== "Response2"

    ```json
    {
      "cases_related_to_created_entities_count": 0,
      "cases_related_to_updated_entities_count": 0,
      "code": 200,
      "created_entity_count": 0,
      "entities": [
        {
          "action": "version",
          "errors": [],
          "id": "32e9fd2c-877a-4700-a06f-bb34e0590ca5",
          "related_cases": [
            {
              "id": "b10c64c2-7fd2-4210-b975-034affb14b57",
              "submitter_id": "TCGA-4G-AAZT"
            }
          ],
          "type": "clinical_supplement",
          "unique_keys": [
            {
              "project_id": "TCGA-CHOL",
              "submitter_id": "nationwidechildrens.org_CHOL.bio.Level_1.428.25.0.tar.gz_nationwidechildrens.org_clinical.TCGA-4G-AAZT-.xml"
            }
          ],
          "valid": true,
          "warnings": []
        }
      ],
      "entity_error_count": 0,
      "message": "Transaction successful.",
      "success": true,
      "transaction_id": 922607,
      "transactional_error_count": 0,
      "transactional_errors": [],
      "updated_entity_count": 0
    }
    ```

### Downloading Files

Files in `file_state = validated` can be downloaded by the submitter using the API or the Data Transfer Tool. This is done in a similar manner as files available in the Data Portal, but will require submission access to the particular project in dbGaP as opposed to downloader access.  File UUIDs can be found in the original upload manifest file, the submission portal, or by API calls.  See [Downloading Files](Downloading_Files.md) for details.

### Deleting Files

Uploaded files must be deleted using a two step process.  First, the file is deleted using the Data Transfer Tool.  See [Deleting Previously Uploaded Data](/Data_Transfer_Tool/Users_Guide/Data_Download_and_Upload/#deleting-previously-uploaded-data) for details.

Second, the file node can be deleted or modified. See [Deleting Entities](#deleting-entities) for details.

## Querying Submitted Data Using GraphQL

### GraphQL Overview

[GraphQL](https://graphql.org/) is a query language that makes it easy to search and retrieve data from graph data structures such as the GDC Data Model.

Unlike the methods outlined in [Search and Retrieval](Search_and_Retrieval.md), which provide access to public releases (or snapshots) of GDC data, the `/graphql` endpoint of GDC Submission API makes it possible for submitters to access "live" data, which provides a real-time view of the state of entities in a project.

>**NOTE:** Access to GDC Submission API GraphQL service is limited to authorized and authenticated submitters. Submitters may only access data in their own project using GraphQL.


### GraphQL IDE

The GDC GraphQL IDE is an instance of [GraphiQL](https://github.com/graphql/graphiql), an in-browser GraphQL IDE that facilitates construction and execution of GraphQL queries. The GDC GraphQL IDE provides tab-completion and syntax checking using schema from the GDC Data Dictionary. It can be found at [https://portal.gdc.cancer.gov/submission/graphiql](https://portal.gdc.cancer.gov/submission/graphiql).

Before interacting directly with the GDC Submission API's GraphQL endpoint, users are encouraged to become familiar with executing queries using the GDC GraphQL IDE.

### GraphQL Endpoint

GDC data submitters can access the GDC Submission API GraphQL endpoint at:

	https://api.gdc.cancer.gov/[API_version/]submission/graphql

where __[API_version/]__ is the optional API version component (see [Getting Started](Getting_Started.md)).

>**NOTE:** An authentication token is required for all requests to the `graphql` endpoint. Queries are restricted to those projects for which the submitter has obtained authorization.


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

In its simplest form, a GraphQL query is a **selection set** (curly brackets) that encloses a set of **fields**. The selection set defines the set of information that is to be retrieved. Furthermore, in GraphQL fields are conceptually equivalent to functions that retrieve additional fields and, in some cases, can take arguments. So each field in a selection set can have its own selection set, thereby creating a nested query structure that can navigate complex data relationships. See [GraphQL Specification](https://spec.graphql.org/) for further details.

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

1. [Escape](https://www.freeformatter.com/json-escape.html) the query using JSON string rules
2. Wrap the query in a ["query" JSON object](#constructing-a-query).
3. Pass the query to the `graphql` endpoint in an HTTP POST request.

Using the `case` and `_case_count` example above as the starting point, the results are as follows:

=== "bare_GraphQL"

    ```GraphQL
    {
    	case (project_id: "TCGA-ALCH", first: 0) {
    		id
    		submitter_id

    	}
    	_case_count (project_id: "TCGA-ALCH")
    }
    ```

=== "escaped_GraphQL"

    ```GraphQL
    {\n\tcase (project_id: \"TCGA-ALCH\", first: 0) {\n\t\tid\n\t\tsubmitter_id\n\n\t}\n\t_case_count (project_id: \"TCGA-ALCH\")\n}
    ```

=== "Query_json"

    ```json
    {
    	"query": "{\n\tcase (project_id: \"TCGA-ALCH\", first: 0) {\n\t\tid\n\t\tsubmitter_id\n\n\t}\n\t_case_count (project_id: \"TCGA-ALCH\")\n}",
    	"variables": null
    }
    ```

=== "Shell_command"

    ```shell
    token=$(<gdc-token-text-file.txt)

    curl --request POST --header "X-Auth-Token: $token" 'https://api.gdc.cancer.gov/v0/submission/graphql' --data-binary @Query_json
    ```

=== "API_Response"

    ```json
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

#### Example: File UUID

GraphQL query to find the file UUID based on file `submitter_id`:

=== "bare_GraphQL"

    ```GraphQL
    {

      submitted_unaligned_reads (project_id: "GDC-INTERNAL", submitter_id: "Blood-00001-aliquot_lane1_barcode23.fastq") {
        id
        submitter_id
        file_name
        project_id
    }
    }
    ```

=== "escaped_GraphQL"

    ```GraphQL
    {
        "query": "{\n \n  submitted_unaligned_reads (project_id: \"GDC-INTERNAL\", submitter_id: \"Blood-00001-aliquot_lane1_barcode23.fastq\") {\n    id\n    submitter_id\n    file_name\n    project_id\n}\n}",
        "variables": null
    }
    ```

=== "Shell"

    ```Shell
    curl --request POST --header "X-Auth-Token: $token" 'https://api.gdc.cancer.gov/v0/submission/graphql' --data-binary @escaped_GraphQL
    ```

=== "Response"

    ```json
    {
      "data": {
        "submitted_unaligned_reads": [
          {
            "file_name": "dummy.fastq",
            "id": "616eab2f-791a-4641-8cd6-ee195a10a201",
            "project_id": "GDC-INTERNAL",
            "submitter_id": "Blood-00001-aliquot_lane1_barcode23.fastq"
          }
        ]
      }
    }
    ```



#### Example: Case Without Diagnosis

GraphQL query for any one case in 'TCGA-LUAD' without Diagnosis information:

=== "bare_GraphQL"

    ```GraphQL
    {
      case (project_id: "TCGA-LUAD", without_links: ["diagnoses"], first: 1) {
          submitter_id
      }
    }
    ```

=== "Response"

    ```json
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

#### Example: Number of Cases Without Diagnosis

GraphQL query for the number of cases in 'TCGA-LUAD' without Diagnosis information:

=== "bare_GraphQL"

    ```GraphQL
    {
      _case_count (project_id: "TCGA-LUAD", without_links: ["diagnoses"])
    }
    ```

=== "Response"

    ```json
    {
      "data": {
        "_case_count": 5
      }
    }
    ```

#### Example: Aliquot State

Query for the `state` of aliquots belonging to case with `submitter_id: "TCGA-ALCH-000001"`:

=== "bare_GraphQL"

    ```GraphQL
    {
      aliquot(with_path_to: {type: "case", submitter_id:"TCGA-ALCH-000001"}) {
        id release_state
      }
    }
    ```

=== "Response"

    ```json
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

#### Example: Aliases

GraphQL query that uses a GraphQL fragment to get specific properties from two portions and give them aliases in the response:

=== "bare_GraphQL"

    ```GraphQL
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

=== "Response"

    ```json
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

#### Example: Biospecimen Tree

GraphQL Query for a case in "TCGA-LUAD" and return a biospecimen tree:

=== "bare_GraphQL"

    ```graphQL
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

=== "Response"

    ```json
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
