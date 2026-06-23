# Appendix C: Format of Submission Requests and Responses

### Format of Submission Request

The general format of JSON objects submitted to the GDC API is as follows:

	{
	    "type": string,
	    "id": string,
	    "submitter_id": string,
	    "<properties>": any type,
	    "<relationship_name>": [
	        {
	            "id": string,
	            "submitter_id": string
	        },
	        ...
	    ]
	}

The request must specify either an `id` or a `submitter_id`.

**`id`**: A string specifying the `id` of the node that the user is creating or updating. This is the persistent GDC UUID4 for the node. If it is preferred to refer to the node using a custom id, users can do so with the `submitter_id` field (described below).

**`submitter_id`**: A string specifying the custom id of the object the user is creating or updating. This is not the official GDC ID for the node.

**`<properties>`**: These key-value pairs will be applied to properties on the referenced node.

**`<relationship_name>`**: A JSON object that specifies a relationship (link) between the node and other nodes. Links are typically established using the `submitter_id` or `id` of the neighboring node.

### Format of API Response to a Submission Request

The following fields are included in all API responses to submission requests.

	{
	  "cases_related_to_created_entities_count": int,
	  "cases_related_to_updated_entities_count": int,
	  "code": int,
	  "created_entity_count": int,
	  "entities": [entities],
	  "entity_error_count": int,
	  "message": string,
	  "success": boolean,
	  "transaction_id": string,
	  "transactional_error_count": int,
	  "transactional_errors": [transactional_errors],
	  "updated_entity_count": int
	}

**`cases_related_to_created_entities_count`**: Number of cases related to the created entities.

**`cases_related_to_updated_entities_count`**: Number of cases related to the updated entities.

**`code`**: The HTTP status code of the response message.

**`created_entity_count`**: Number of entities created.

**`entities`**: A list of entities of the form:

	{
	  "action": string,
	  "errors": [entity_errors],
	  "id": string,
	  "related_cases": [object],
	  "type": string,
	  "unique_keys": [unique_keys],
	  "valid": boolean,
	  "warnings": [object]
	}

*`entity_errors`*: A list of errors that occurred while parsing, validating, or performing a CRUD operation on a
specific entity. Entity errors are of the form:

	{
		"keys": [string],
		"message": string
	}


*`unique_keys`*: Properties, or combinations of properties, that can be used to uniquely identify the node in the GDC.  Unique_keys are of the form:


	{
		"project_id": string,
		"submitter_id": string
	}


**`entity_error_count`**: Number of entities that were not successful.

**`message`**: A human-readable message describing the transaction.

**`success`**: A boolean value stating whether the transaction was successful. If the value is False, then no changes will be made to the database.

**`transaction_id`**: A string specifying the transaction id.

**`transactional_error_count`**: Number of transactional errors that occurred.

**`transactional_errors`**: A list of transactional errors that have occurred. These errors are errors that are not specific to an individual entity. Transactional errors are of the form:

	{
		"message": string
	}

**`updated_entity_count`**: Number of existing entities updated by the transaction.


### Error Types

**`EntityNotFoundError`** A referenced entity was not found among existing entities and entities specified in the transaction.

**`MissingPropertyError`** A required property was not provided.

**`ValidationError`** A provided property did not pass a validation test.

## Status Messages

API responses will contain a status for each entity specified in the request:

**`success`**: The desired transaction was sucessful and the entity's state was modified in the GDC.

**`valid`**: The desired transaction was not sucessful, but the trasaction was not aborted because of this entity.

**`error`**: The desired transaction was not sucessful, and the transaction was aborted because of this entity. This entity did not pass validation or an internal error occured when attempting to complete the transaction. The error state will be accompanied by a list of errors recorded about the entity (see label-error-messages).

>**Note:** GDC API requests are transactional. An error with processing a node specified in the transaction will abort the transaction and will result in no changes being applied for any node involved in the transaction.
