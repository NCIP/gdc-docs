# BAM Slicing

The GDC API provides remote BAM slicing functionality that enables downloading of specific parts of a BAM file instead of the whole file. This functionality can be accessed at the `slicing` endpoint, using a syntax similar to that of widely used bioinformatics tools such as `samtools`.

## About the slicing endpoint

The `slicing` endpoint accepts HTTP GET requests in the form of a URL, and HTTP POST requests that carry a JSON payload. POST requests are more appropriate in cases where query parameters make the GET URL very long.

The response will be a BAM-formatted file containing the header of the source BAM file, as well as any alignment records that are found to overlap the specified regions, sorted by chromosomal coordinate.

Please note the following:

* The functionality of this API differs from the usual functionality of `samtools` in that alignment records that overlap multiple regions will not be returned multiple times.
* A request with no region or gene specified will return the BAM header, which makes it easy to inspect the references to which the alignment records were aligned.
* A request for regions that are not included in the source BAM is not considered an error, and is treated the same as if no records existed for the region.
* Examples provided for BAM slicing functionality are intended for use with GDC harmonized data (i.e. BAM files available in the GDC Data Portal). Slicing of unharmonized BAM files (i.e. BAM files in the GDC Legacy Portal) is not supported.
* Bam slicing does not create an associated bam index (.bai) file.  For applications requiring a .bai file users will need to generate this file from the bam slice using a tool and command such as `samtools index`.

### Query Parameters

The following query parameters and JSON fields are supported:

| Description | Query Parameter | JSON Field | Query format |
|---|---|---|---|
| entire chromosome, or a position or region on the chromosome, specified using chromosomal coordinates | region | regions | region=<chr>(:<start>(-<stop>)?)?</stop></start></chr> |
| region specified using a [HGNC](http://www.genenames.org/) / [GENCODE v22](http://www.gencodegenes.org/) gene name |  gencode | gencode | gencode=<gene_name> |

>**NOTE:** The successfully sliced BAM will contain all reads that overlap (entirely or partially) with the specified region or gene. It is possible to specify an open-ended region, e.g. `chr2:10000`, which would return all reads that (completely or partially) overlap with the region of chromosome 2 from position 10,000 to the end of the chromosome.

### JSON Schema

JSON payloads can be syntactically verified using the following JSON schema:

	{
	  "$schema": "http://json-schema.org/schema#",
	  "type": "object",
	  "properties": {
	    "regions": {
	      "type": "array",
	      "items": {
	        "type": "string",
	        "pattern": "^[a-zA-Z0-9]+(:([0-9]+)?(-[0-9]+)?)?$"
	      }
	    },
	    "gencode": {
	      "type": "array",
	      "items": {
	        "type": "string"
	      }
	    }
	  }
	}



## Examples: Specifying a region

The following two requests are examples of BAM slicing using region(s).

```Regions_GET

token=$(<gdc-token-text-file.txt)

curl --header "X-Auth-Token: $token" 'https://api.gdc.cancer.gov/slicing/view/df80679e-c4d3-487b-934c-fcc782e5d46e?region=chr1&region=chr2:10000&region=chr3:10000-20000' --output get_regions_slice.bam
```
```Regions_POST_Payload
{
    "regions": [
        "chr1",
        "chr2:10000",
        "chr3:10000-20000"
    ]
}
```
```Regions_POST
token=$(<gdc-token-text-file.txt)

curl --header "X-Auth-Token: $token" --request POST https://api.gdc.cancer.gov/slicing/view/9ca90dfa-e62f-4f9c-9946-dfcecfd3ca4d --header "Content-Type: application/json" -d@Payload --output post_regions_slice.bam
```
```Response
Response:
HTTP/1.1 206

<bam_data_stream>
```

## Examples: Specifying a gene

The following two requests are examples of BAM slicing using HGNC / GENCODE v22 gene name(s).

```Gencode_GET
token=$(<gdc-token-text-file.txt)

curl --header "X-Auth-Token: $token" 'https://api.gdc.cancer.gov/slicing/view/df80679e-c4d3-487b-934c-fcc782e5d46e?gencode=BRCA1' --output get_brca1_slice.bam
```
```Gencode_POST_Payload
{
    "gencode": [
        "BRCA1",
        "BRCA2"
    ]
}
```
```Gencode_POST
curl --header "X-Auth-Token: $token" --request POST https://api.gdc.cancer.gov/slicing/view/df80679e-c4d3-487b-934c-fcc782e5d46e --header "Content-Type: application/json" -d@Payload --output post_brca12_slice.bam
```
```Response
Response:
HTTP/1.1 206

<bam_data_stream>
```
## Examples: Specifying unmapped reads

Unmapped reads are found in GDC BAM files. You may request these reads by using the following commands.

```GET
token=$(<gdc-token-text-file.txt)

curl --header "X-Auth-Token: $token" 'https://api.gdc.cancer.gov/slicing/view/dc87e1b8-d8b7-4837-88ea-fb7f017b3c69?region=unmapped' --output get_regions_slice.bam
```
```POST_Payload
{
    "regions": [
        "unmapped"
    ]
}
```
```POST
curl --header "X-Auth-Token: $token" --request POST https://api.gdc.cancer.gov/slicing/view/dc87e1b8-d8b7-4837-88ea-fb7f017b3c69 --header "Content-Type: application/json" -d@Payload --output get_regions_slice.bam
```
```Response
Response:
HTTP/1.1 206

<bam_data_stream>
```



After downloading, the sliced BAM file can be converted to SAM using the following command if `samtools` is installed on the user's system:

	samtools view -h brca1_slice.bam -o brca1_slice.sam

## Errors

When slicing cannot be performed, the GDC API will provide JSON error responses and HTTP error codes.

### JSON Error Responses

JSON error responses have the following structure:

	{
	    "error": "<error-message>"
	}

For example, when making a request for a protected BAM without supplying a GDC authentication token:

```Shell
curl https://api.gdc.cancer.gov/v0/slicing/view/15b0bf8e-ff20-41ab-8366-a495c11b30be
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

Early truncation can arise when connection is interrupted or when slicing fails due to BAM corruption.
