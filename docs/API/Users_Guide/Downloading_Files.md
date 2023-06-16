# Downloading Files

The GDC API implements file download functionality using `data` and `manifest` endpoints. The `data` endpoint allows users to download files stored in the GDC by specifying file UUID(s). The `manifest` endpoint generates a download manifest file that can be used with the GDC Data Transfer Tool to transfer large volumes of data.

>**Note:** Downloading controlled access data requires the use of an authentication token. See [Getting Started: Authentication](Getting_Started.md#authentication) for details.

## Data endpoint

To download a file, users can pass UUID(s) to the `data` endpoint.  If a single UUID is provided, the API will return the associated file. If a comma-separated list of UUIDs is provided, the API will return an archive file containing the requested files.

The `data` endpoint supports GET and POST requests as demonstrated in the following examples.


### Downloading a Single File using GET

This example demonstrates downloading a single file from the GDC. Here we pass the file's UUID to the `data` endpoint with a GET request.

```shell
curl --remote-name --remote-header-name 'https://api.gdc.cancer.gov/data/5b2974ad-f932-499b-90a3-93577a9f0573'
```
```Output
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100 6111k  100 6111k    0     0   414k      0  0:00:14  0:00:14 --:--:--  412k

```
### Related Files

If the `related_files=true` parameter is specified, the following related files, if available, will be included in the download package by the GDC API:

* BAM index files (BAI files)
* VCF index files (TBI files)

For example, this request will download a BAM file and its associated BAI file:

```shell
curl --remote-name --remote-header-name -H "x-auth-token: $token" "https://api.gdc.cancer.gov/data/f587ef82-acbe-44f9-ad5a-6207e148f61f?related_files=true"
```
```Output
% Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                               Dload  Upload   Total   Spent    Left  Speed
100 63.4M    0 63.4M    0     0  7541k      0 --:--:--  0:00:08 --:--:--  9.9M

```


### Downloading Multiple Files using GET

This example demonstrates downloading multiple files from the GDC using a GET request. The GDC API returns a `.tar.gz` archive containing the downloaded files.

```shell
curl --remote-name --remote-header-name 'https://api.gdc.cancer.gov/data/e3228020-1c54-4521-9182-1ea14c5dc0f7,18e1e38e-0f0a-4a0e-918f-08e6201ea140'
```
```Output
% Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                               Dload  Upload   Total   Spent    Left  Speed
100  287k    0  287k    0     0  30131      0 --:--:--  0:00:09 --:--:-- 42759

```

>**Note:** This method supports downloading a limited number of files at one time. To download a large number of files, please use [POST](#downloading-multiple-files-using-post).

#### Downloading an Uncompressed Group of Files

If the `?tarfile` parameter is specified to a data endpoint download query all files requested in the download string will be bundled in a single tar file rather than a tar.gz file which is the default behavior.  
```shell
curl --remote-name --remote-header-name 'https://api.gdc.cancer.gov/data/1da7105a-f0ff-479d-9f82-6c1d94456c91,77e73cc4-ff31-449e-8e3c-7ae5ce57838c?tarfile'
```      

### Downloading Multiple Files using POST

The following two examples demonstrate downloading multiple files from the GDC using a POST request that contains a payload in one of two formats: percent-encoded form data or JSON. The GDC API returns a `.tar.gz` archive containing the downloaded files.

#### POST request with form data payload

POST requests that carry a payload of percent-encoded form data must include the HTTP header `Content-Type: application/x-www-form-urlencoded`.

The payload is a string in the following format:

	ids=UUID1&ids=UUID2&ids=UUID3...

where UUID# corresponds to the UUIDs of the files to be downloaded.

In this example we use `curl` to download a set of files from the GDC Data Portal. The payload is stored in a plain text file named `Payload`; `curl` includes the `Content-Type: application/x-www-form-urlencoded` header by default.

```Payload
ids=59eb3fc5-9172-4828-8dec-0d9988073103&ids=869b7d7c-ff35-482a-aa8d-1a8675c161d3&ids=b8ffff40-aa0e-4534-b05f-9311f16c2f6b&ids=51e14969-30a7-42d9-8168-4a5ea422ca4a&ids=adcfc856-990b-40fc-8f1e-67dfc2343fb7&ids=7f1e9aee-eb4e-4c79-8626-b603c9be124d&ids=62a8feb5-c660-4261-bcd6-67fbb79bb422
```
```Shell
curl --remote-name --remote-header-name --request POST 'https://api.gdc.cancer.gov/data' --data @Payload
```
```Output
% Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                Dload  Upload   Total   Spent    Left  Speed
100 6804k    0 6804k  100   286   245k     10  0:00:28  0:00:27  0:00:01  357k
```

#### POST request with JSON payload

POST requests that carry a JSON payload must include the HTTP header `Content-Type: application/json`.

The payload is a string in the following format:

	{
	    "ids":[
	        "UUID1",
	        "UUID2",
					...
	        "UUID3"
	    ]
	}


where UUID# corresponds to the UUIDs of the files to be downloaded.

In this example we use `curl` to download a set of files from the GDC Portal; the payload is stored in a plain text file named `Payload`.


```Payload
{
    "ids":[
        "0451fc55-33ef-4151-a68c-cac59be716dc",
        "0cc3d450-2c60-4cb0-a073-d92dc979fa5e",
        "0de9bc40-3ef8-4fe7-b7d6-80a9339b0bf8",
        "0f8d8202-a1ca-4ea1-98b2-c20a6b08479a"
    ]
}
```
```Shell
curl --remote-name --remote-header-name --request POST --header 'Content-Type: application/json' --data @request.txt 'https://api.gdc.cancer.gov/data'
```
```Output
% Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                Dload  Upload   Total   Spent    Left  Speed
100 5878k    0 5878k  100   205   290k     10  0:00:20  0:00:20 --:--:--  198k
```

### Downloading Controlled-access Files

To download controlled-access files, a valid authentication token must be passed to the GDC API using the `X-Auth-Token` HTTP header:

```shell
token=$(<gdc-token-text-file.txt)

curl --remote-name --remote-header-name --header "X-Auth-Token: $token" 'https://api.gdc.cancer.gov/data/be6d269d-4305-4643-b98e-af703a067761'
```
```Output
% Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                             Dload  Upload   Total   Spent    Left  Speed
100 65.8M  100 65.8M    0     0   271k      0  0:04:08  0:04:08 --:--:--  288k
```

## Manifest endpoint

The `manifest` endpoint generates a download manifest file that can be used with the [GDC Data Transfer Tool](../../Data_Transfer_Tool/Users_Guide/Getting_Started.md). The Data Transfer Tool is recommended for transferring large volumes of data. The GDC API can also generate a download manifest from a list of results that match a [Search and Retrieval](Search_and_Retrieval.md) query. To do this, append `&return_type=manifest` to the end of the query. Note that the "size" parameter does not work with the manifest endpoint and will return the entire set of files.


### Using the manifest endpoint

The `manifest` endpoint allows users to create a download manifest, which can be used with the GDC Data Transfer Tool to download a large volume of data. The `manifest` endpoint generates a manifest file from a comma-separated list of UUIDs.

```shell
curl --remote-name --remote-header-name 'https://api.gdc.cancer.gov/v0/manifest/a751cc7e-d2ff-4e9a-8645-09bf12612f1a,9c97e3fe-1610-4a92-9a24-ab3b9e4000e2'
```
```Output
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100   274  100   274    0     0   1042      0 --:--:-- --:--:-- --:--:--  1041

```

The `manifest` endpoint also supports HTTP POST requests in the same format as the `data` endpoint; see [above](#post-request-with-json-payload) for details.


### Using return_type=manifest

Alternatively, users can create a manifest by appending `&return_type=manifest` to a [Search and Retrieval](Search_and_Retrieval.md) query. In this example, we generate a download manifest for RNA-seq data files from solid tissue normal samples, that are part of the TCGA-KIRC project:

```Shell
curl --remote-name --remote-header-name 'https://api.gdc.cancer.gov/files?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B%22field%22%3A%22experimental_strategy%22%2C%22value%22%3A%5B%22RNA-Seq%22%5D%7D%7D%2C%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%5B%22TCGA-KIRC%22%5D%7D%7D%2C%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B%22field%22%3A%22cases.samples.sample_type%22%2C%22value%22%3A%5B%22Solid+Tissue+Normal%22%5D%7D%7D%5D%7D&return_type=manifest'
```
```Output
% Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                Dload  Upload   Total   Spent    Left  Speed
100  100k    0  100k    0     0   277k      0 --:--:-- --:--:-- --:--:--  282k
```
