# Submission

## Overview

Using the following methods, users submitting to the GDC can create, delete and update entities and relationships in the GDC data model.

## GDC Dictionary

Requests to the submission API must adhere to the schemas defined in the [GDC Data Dictionary](https://www.github.com/NCI-GDC/gdcdictionary).

## Working with Entities

### Query Format

When updating, creating, or deleting entities in the GDC, users need to specify the entity type, the entity id, any relationships the entity has to parent entities from which it was derived, and any properties (required and optional as defined by the entity schema). The structure for each entity should look as follows:

```json
{
    "type": string,
    "id": string,
    "submitter_id": string,
    "<entity_property_keys>": any type,
    "<relationship_name>": [
        {
            "id": string,
            "submitter_id": string
        },
        ...
    ]
}
```

The request must specify either an `id` or a `submitter_id`.

**`id`** : A string specifying the id of the entity the user is creating, updating, or deleting. This is the official GDC UUID for the entity. If it is prefered to refer to the entity using a custom id, users can do so with the `submitter_id` field (described below).

**`submitter_id`** : A string specifying the custom id of the object the user is creating, updating or deleting. This is not the official GDC ID for the entity.

**`<entity_property_keys>`** : All keys except for `id` and `submitter_id` will be treated as properties keys. These key value pairs will be used as properties on referenced entity.

**`<relationship_name>`** : A JSON object specifying a relationship. The value for this is a JSON object specifying either the submitter_id or the id of the neighboring entity.

### Response Format

The following fields are included in all API responses to submission requests.

```json
{
	"code": int,
	"created_entity_count": int,
	"entities": [object],
	"entity_error_count": int,
	"message": string,
	"success": boolean,
	"transactional_error_count": int,
	"transactional_errors": [transactional_error],
	"updated_entity_count": int
}
```

**success** A boolean value stating whether the transaction was successful. If the value is False, then no changes will be made to the database.

**code** The HTTP status code of the response message. A human readable summary of the transaction results.

**transactional_errors** A list of transactional errors that have occurred. These errors are errors that are not specific to
an individual entity. Transactional errors are of the form:

```json
{
	"message": string
}
```

**`transactional_error_count`** A count of the number of transactional errors that occured.

**`entity_error_count`** A count of the number of entities that were not successful.

**entities** A list of entities of the form:

```json
{
	"submitter_id": string,
	"errors": [entity_errors],
	"id": string,
	"valid": boolean,
	"type": string
}
```

**entity_errors**

A list of errors that occurred while parsing, validating, or performing a CRUD operation on a
specific entity. Entity errors are of the form:

```json
{
	"keys": [string],
	"message": string
}
```

For a listing of the types of errors, see Creating Entities.

**`created_entity_count`** The number of entities created by the transaction.

**`updated_entity_count`** The number of existing entities updated by the transaction.

### Creating Entities

Entities can be created via both the POST and PUT HTTP methods.

**When to use POST:** Submitters should use the POST method
when creating new entities. The POST method will ensure that the entities in the posted transaction have not
previously been created. This helps to avoid inadvertent updates that can occur when using PUT.

**When to use PUT:** Submitters should use the PUT method to update
existing entities and/or create new ones. The PUT method will create entities that have not been previously created and update those that have been previously created. The response body will specify whether each entity in the transaction was created or updated.

Requests to create/update entities are transactional, meaning that if a single entity is invalid, the transaction will fail and no changes will be made to the database.

`POST /v0/submission/<program>/<project>/` This endpoint is used to create GDC entities. Using the POST on a project’s endpoint will
create any valid entities specified in the request body.

**Parameters**

* program (str) – The program to which the submitter belongs and the context in
which the API request is valid. The program is the human-readable name, e.g.
TCGA.

* project (str) – The project to which the submitter belongs and the context in
which The API request is valid. The project is the human-readable code, e.g.
BRCA.

**Example Usage:**

The following example will:

1. Create a new aliquot *aliquot-1*.

2. Specify that *aliquot-1* was derived from analyte *analyte-1*.

3. Specify that *analyte-1* was derived from existing portion *portion-1*.

```
	POST /v0/submission/program1/project1/ HTTP/1.1
	Host: example.com
	Content-Type: application/json
	X-Auth-Header: MIIDKgYJKoZIhvcNAQcC...
	Accept: application/json
```

```json
[
    {
        "type": "analyte",
        "portions": {
            "submitter_id": "portion-1"
        },
        "well_number": null,
        "analyte_type": "DNA",
        "submitter_id": "analyte-1",
        "amount": 10.98,
        "a260_a280_ratio": null,
        "concentration": 0.14,
        "spectrophotometer_method": "PicoGreen",
        "analyte_type_id": "D"
    }, {
        "type": "aliquot",
        "analytes": {
            "submitter_id": "analyte-1"
        },
        "submitter_id": "aliquot-1",
        "amount": null,
        "source_center": "23",
        "concentration": 0.07
    }
]
```

**Example Successful Response:**

```
HTTP/1.1 201 CREATED
Content-Type: application/json
```

```json
{
    "code": 201,
    "created_entity_count": 1,
    "entities": [
        {
            "submitter_id": "analyte-1",
            "errors": [],
            "id": "2e1429d5-b2ec-4c02-93ac-207d10b1193c",
            "valid": true,
            "type": "analyte"
        }, {
            "submitter_id": "aliquot-1",
            "errors": [],
            "id": "6a30b20a-1e38-4c16-8c16-c03ab30f7a11",
            "valid": true,
            "type": "aliquot"
        }
    ],
    "entity_error_count": 0,
    "message": "Transaction successful.",
    "success": true,
    "transactional_error_count": 0,
    "transactional_errors": [],
    "updated_entity_count": 0
}
```

**Successful Example Response:**
In the successful example response, the analyte will be created, assigned ID 2e1429d5... and
linked to portion *portion-1*. The aliquot will be created, assigned an ID and linked to analyte
*analyte-1*. Note that the portion *portion-1* referenced by *analyte-1* above is not included
in the transaction. Part of the validation for the creation of any entity is to check if:

1. The entity it was derived from exists in the current transaction. If the parent entity was in the
transaction, verify that any information provided does not conflict with the existing version.

2. If the parent entity was not in the transaction, verify that it already exists in the system.

Just as the portion referenced by *portion-1* was submitted in a previous transaction. This example could have
been split into two transactions, the first creating the aliquot and the second creating the file.

**Note:** GDC will not allow submitters to create entities that do not have relationships to either existing entities in GDC or entities included in the submission transaction. For example, users cannot submit an aliquot if users are not submitting/have not previously
submitted the sample, portion, or analyte from which it was derived. This rule
applies to deleting entries as well (7.2.6 - Deleting Entities.) Creation, updates and deletions of Programs and Projects are an administrative functions handled by GDC.

**Example Bad Request:**

```
POST /v0/submission/program1/project1 HTTP/1.1
Host: example.com
Content-Type: application/json
X-Auth-Header: MIIDKgYJKoZIhvcNAQcC...
Accept: application/json
```

```json
[
    {
        "type": "analytes",
        "portions": {
            "submitter_id": "portion-1"
        },
        "well_number": null,
        "analyte_type": "DNA",
        "submitter_id": "analyte-1",
        "amount": 10.98,
        "a260_a280_ratio": null,
        "concentration": 0.14,
        "spectrophotometer_method": "PicoGreen",
        "analyte_type_id": "D"
    }, {
        "type": "aliquot",
        "analytes": {
            "submitter_id": "analyte-1"
        },
        "submitter_id": "aliquot-1",
        "amount": null,
        "source_center": "23",
        "concentration": 0.07
    }
]
```

**Example Error Result:**

```
HTTP/1.1 400 BAD REQUEST
Content-Type: application/json
```

```json
{
    "code": 400,
    "created_entity_count": 0,
    "entities": [
        {
            "submitter_id": "analyte-1",
            "errors": [
                {
                    "keys": ["type"],
                    "message": "Invalid entity type: aliquots. Did you mean 'aliquot'?"
                }
            ],
            "id": "2e1429d5-b2ec-4c02-93ac-207d10b1193c",
            "valid": false,
            "type": "analyte"
        }, {
            "submitter_id": "aliquot-1",
            "errors": [],
            "id": "6a30b20a-1e38-4c16-8c16-c03ab30f7a11",
            "valid": true,
            "type": "aliquot"
        }
    ],
    "entity_error_count": 1,
    "message": "Transaction aborted due to 1 invalid entity.",
    "success": false,
    "transactional_error_count": 0,
    "transactional_errors": [],
    "updated_entity_count": 0
}
```
**Unsuccessful Example Response:**

In the second example response, the API returned error code 400 and each entity with a list of errors.

The GDC API will also return a list of all errors by entity.

### Retrieving Entities

`GET /v0/submission/<program>/<project>/entities/entity_id_string.` This endpoint is for retrieving existing GDC entities by ID. For more advanced querying on entities or retrieving set of entities, the GraphQL endpoint described in Section 6.4 is recommended.

The return type of a GET on this endpoint is a JSON array containing JSON object elements, each
corresponding to a provided ID. Return results are unordered.

If any ID is not found in the database, a status code of 404 is returned with the missing IDs.

**Parameters**

* **program** (str) – The program to which the case belongs and the context in
which the API request is valid. The program is the human-readable name, e.g.
TCGA.

* **project** (str) – The project to which the case belongs and the context in
which The API request is valid. The project is the human-readable code, e.g.
BRCA.

* **ids** (str) – A comma separated list of ids specifying the entities to retrieve. These
ids may be official GDC ids or project unique submitter_id.

### Updating Entities

`PUT /v0/submission/<program>/<project>/` This endpoint is used to update/create GDC entities. Using the PUT method on a project’s endpoint
will, for any valid entities specified in the request body, create those that do not exist and update those
that do.

**Parameters**

* program (str) – The program to which the case belongs and the context in
which the API request is valid. The program is the human-readable name (e.g.
TCGA).

* project (str) – The project to which the case belongs and the context in
which the API request is valid. The project is the human-readable code (e.g.
BRCA).

The request body syntax is the same as the POST method for the same endpoint.

### Deleting Entities


`DELETE /v0/submission/<program>/<project>/entities/ids`.

The above endpoint is used to delete existing GDC entities. Using DELETE on the a project’s endpoint will completely delete an entity.

The GDC does not allow deletions or creations that would leave nodes without parents (i.e. nodes that do not have an entity from which they were derived). To prevent catastrophic mistakes the automatic cascading of deletes is not allowed.

To inform the user which entities must be deleted for the target entity to be deleted, the API will respond with a list of entities that must be deleted prior to deleting the target entity.

**Parameters**

* **program** (str) – The program to which the submitter belongs and the context in
which the API request is valid. The program is the human-readable name (e.g.
TCGA).

* **project** (str) – The project to which the submitter belongs and the context in
which the API request is valid. The project is the human-readable code (e.g. BRCA).

* **ids** (str) – A comma separated list of ids specifying the entities to delete. These ids must be official GDC ids.

### Error Types

**EntityNotFoundError** A referenced entity was not found. This includes both the transaction and the datamodel.

**MissingPropertyError** A required property was not provided.

**ValidationError** A provided property did not pass a validation test.

### Status Messages

API responses will contain a status for each entity specified in the request:

**Note:** Since GDC API requests are transactional, either all entities will be processed successfully, resulting in a success status message or none will (i.e. one invalid entity will result in the entire transaction being aborted). Information contained in the status code '**error**' can help users resolve issues.

**success**: The desired transaction was sucessful and the entity’s state was modified in the database.

**valid**: The desired transaction was not sucessful, but the trasaction was not aborted because of this entity. Had all other entities in this transaction been valid and there were no internal errors, then the status
of this entity would succeed.

**error**: The desired transaction was not sucessful, and the transaction was in part aborted because of this
entity. This entity did not pass validation or an internal error occured when attempting to complete
the transaction. The error state will be accompanied by a list of errors recorded about the entity
(see label-error-messages).

## GraphQL

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

#### Sample GraphQL query

The following is a GraphQL query for a case in project
TCGA-LAML that returns a JSON document containing the `submitter_id` of
the case and of its samples.


```JavaScript
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
