# BAM Slicing

The GDC API provides BAM slicing functionality via the `slicing` endpoint. It enables remote slicing of BAM files, using a syntax similar to that of widely used bioinformatics tools such as `samtools`.

**NOTE:** Authentication is handled in the same manner as with data download,
via the `X-Auth-Token` header. See [Downloading Controlled-access Files](Downloading_Files.md#downloading-controlled-access-files)

## GET and POST Queries

The `slicing` endpoint accepts HTTP `GET` requests in the form of a URL, and HTTP `POST` requests that carry a JSON payload. `POST` requests are more appropriate in cases where query parameters make the `GET` URL very long.


```GET
GET /v0/slicing/view/fd2386b3-32f5-4aea-a55c-af24f5d1e650?region=chr1&region=chr2:1000&region=chr3:1000-2000 HTTP/1.1
X-Auth-Token: <gdc-user-token>
```
```POST
POST /v0/slicing/view/fd2386b3-32f5-4aea-a55c-af24f5d1e650 HTTP/1.1
Content-Type: application/json
X-Auth-Token: <gdc-user-token>

{
    "regions": [
        "chr1",
        "chr2:1000",
        "chr3:1000-2000"
    ]
}
```
```Response
Response:
HTTP/1.1 206

<bam_data_stream>
```

## Response
The response will be a BAM-formatted file containing the header of the source BAM file, as well as any alignment records that are found to overlap the specified regions, sorted by chromosomal coordinate.

**NOTE:** The functionality of this API differs from the usual functionality of `samtools` in that alignment records that overlap multiple regions will not be returned multiple times.

## Query Parameters
The following query parameters and JSON fields are supported:

| Query Parameter   | JSON Field      |
|----------|:-------------:|
| **region:** a chromosomal-coordinate region | region=<chr>(:<start>(-<stop>)?)?</stop></start></chr> |
| **bai:** an overriding GDC BAI ID to be used for region lookup |    bai=<gdc_bai_id></gdc_bai_id>   |  


JSON payloads can be syntactically verified using the following JSON schema:

	{
	    "type": "object",
	    "properties": {
	        "bai": {
	            "type": "string"
	        },
	        "regions": {
	            "type": "array",
	            "items": {
	                "type": "string",
	                "pattern": "^[a-zA-Z0-9]+(:[0-9]+(-[0-9]+)?)?$"
	            }
	        }
	    },
	    "additionalProperties": false
	}


A request with no regions specified will return the BAM header, which makes it easy to inspect the references to which the alignment records were aligned.

**NOTE:** A request for regions that are not included in the source BAM is not considered an error, and is treated the same as if no records existed for the region.

## Errors

When slicing cannot be performed, the GDC API will provide JSON error responses and HTTP error codes.

### JSON Error Responses

JSON error responses have the following structure:

	{
	    "error": "<error-message>"
	}

For example, when making a request for a protected BAM without supplying a GDC authentication token:

```Shell
curl https://api-gdc.nci.nih.gov/v0/slicing/view/fd2386b3-32f5-4aea-a55c-af24f5d1e650
```
```Response
HTTP/1.1 403 FORBIDDEN
{
    "error": "Please specify a X-Auth-Token"
}
```

### HTTP error codes

Potential HTTP error codes include:

Error Code | Description
---------- | -------
400 | Bad Request -- The regions specified are malformed
403 | Unauthorized -- The user could not be authenticated
403 | Unauthorized -- The user is not authorized for access to the source BAM
404 | Not Found -- No BAM is specified
404 | Not Found -- No BAM can be found for the specified GDC BAM ID
504 | BAI Not Found -- No BAI can be found for the BAM


### Transfer Errors

In the case that an error occurs during transfer of the resulting BAM, the BGZF EOF marker will not be present. This early truncation of the BAM file will cause errors if the file is used as input to other programs. For example, `samtools` will provide the error "EOF marker is absent".

Early truncation can arise when:

- Connection is interrupted
- Slicing fails due to BAM corruption
