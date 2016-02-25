# BAM Slicing
GDC exposes a BAM slicing API that allows for remote slicing of BAM format
objects, utilizing a similar syntax to widely used bioinformatics tools, such
as samtools.
<aside class="notice">**NOTE:** Authorization to the API is handled in the same manner as with data download,
via the ```X-Auth-Token``` header.See <a href=#controlled-data-access>Controlled Data Access</a></aside>
## GET Endpoint
The GDC slicing API is exposed at ```/v0/slicing/view/<gdc-bam-id>```. The
following is an example request to this endpoint and the expected response:

```html
Request:
GET /v0/slicing/view/fd2386b3-32f5-4aea-a55c-af24f5d1e650?region=chr1&region=chr2:1000&region=chr3:1000-2000 HTTP/1.1
X-Auth-Token: <gdc-user-token>

Response:
HTTP/1.1 206

<bam_data_stream>
```
## POST Endpoint
In cases where query parameters become unwieldy, users may use a secondary endpoint that is
available that accepts JSON. The following is an identical request made but
with a JSON payload:

```html
Request:
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

Response:
HTTP/1.1 206

<bam_data_stream>
```
## Response
The response will be a BAM-formatted file containing the header of the source
BAM file, as well as any alignment records that are found to overlap the
specified regions, in chromosomal-coordinate sorted order.

<aside class="notice">**NOTE:** The functionality of this API differs from the usual functionality of
samtools in that alignment records that overlap multiple regions will not be
returned multiple times.</aside>
## Query Parameters
The following query parameters and JSON fields are supported:

| Query Parameter   | JSON Field      |
|----------|:-------------:|
| **region:** a chromosomal-coordinate region | region=<chr>(:<start>(-<stop>)?)?</stop></start></chr> | 
| **bai:** an overriding GDC BAI ID to be used for region lookup |    bai=<gdc_bai_id></gdc_bai_id>   |  

>JSON payloads can be syntactically verified using the following JSON Schema:

```json
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
```

A request with no regions specified will simply return the BAM header, which
can be used by the client to quickly inspect the references to which the alignment
records were aligned.

<aside class="notice">**NOTE:** A request for regions that are not included in the source BAM is not
considered an error, and is treated the same as if no records existed for
the region.</aside>
## Errors
>The JSON responses for errors have the following structure:

```json
{
    "error": "<error-message>"
}
```
>The response will contain error messages produced from
the data access endpoints for the GDC API. For example, when making a request
for a protected BAM without supplying a GDC authentication token:

```
curl https://api-gdc.nci.nih.gov/v0/slicing/view/fd2386b3-32f5-4aea-a55c-af24f5d1e650

HTTP/1.1 403 FORBIDDEN
{
    "error": "Please specify a X-Auth-Token"
}
```
The API will provide errors in the form of error codes and JSON when slicing
can not be performed. Such potential errors include:

Error Code | Description
---------- | -------
400 | Bad Request -- The regions specified are malformed
403 | Unauthorized -- The user could not be authenticated 
403 | Unauthorized -- The user is not authorized for access to the source BAM
404 | Not Found -- No BAM is specified 
404 | Not Found -- No BAM can be found for the specified GDC BAM ID
504 | BAI Not Found -- No BAI can be found for the BAM



In the case that an error occurs during transfer of the resulting BAM, the
BGZF EOF marker will not be present. This early truncation of the BAM will cause errors if used as input to other programs. For example, samtools will provide the error "EOF marker is absent". Early truncation can arise when:

- Connection is interrupted
- Slicing fails due to BAM corruption
