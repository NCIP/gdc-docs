# Downloading Files

The GDC API implements file download functionality using `data` and `manifest` endpoints. The `data` endpoint allows users to download files stored in the GDC by specifying file UUID(s). The `manifest` endpoint generates a download manifest file that can be used with the GDC Data Transfer Tool to transfer large volumes of data.

>**Note:** Downloading controlled access data requires the use of an authentication token. See [Getting Started: Authentication](Getting_Started.md#authentication) for details.

>**Note:** Requests to download data from the GDC Legacy Archive may be directed to `legacy/data` or `data`. See [Getting Started: Legacy Archive](Getting_Started.md#gdc-legacy-archive) for details.

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
curl: Saved to filename '14-3-3_beta-R-V_GBL1112940.tif'
```
### Related Files

If the `related_files=true` parameter is specified, the following related files, if available, will be included in the download package by the GDC API:

* BAM index files (BAI files)
* Metadata files (such as SRA XML or MAGE-TAB files)

For example, this request will download a legacy copy number segmentation file and its associated MAGE-TAB metadata file:

```shell
curl --remote-name --remote-header-name 'https://api.gdc.cancer.gov/data/7efc039a-fde3-4bc1-9433-2fc6b5e3ffa5?related_files=true'
```
```Output
% Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
															 Dload  Upload   Total   Spent    Left  Speed
100 65353    0 65353    0     0  65353      0 --:--:-- --:--:-- --:--:--  102k
curl: Saved to filename 'gdc_download_20180830_131817.826097.tar.gz'
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
curl: Saved to filename 'gdc_download_064d1aa8cc8cbab33e93979bebbf7d6af2d6a802.tar.gz'
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

In this example we use `curl` to download a set of files from the GDC Legacy Archive. The payload is stored in a plain text file named `Payload`; `curl` includes the `Content-Type: application/x-www-form-urlencoded` header by default.

```Payload
ids=556e5e3f-0ab9-4b6c-aa62-c42f6a6cf20c&ids=e0de63e2-02f3-4309-9b24-69f4c24e85fc&ids=f1a06178-2ec2-4b06-83f3-3aedac332cfe&ids=11a8aca0-c8e6-4ff8-8ab6-fe18a1b8ba82&ids=69a69c84-00de-45ff-b397-fd2b6713ed4f&ids=9ec48233-395d-401e-b205-951c971f8dd4&ids=93129547-378c-4b69-b858-532abfff678e&ids=8d4277e9-a472-4590-886d-24dc2538ea65&ids=6733b412-56da-4f1c-a12b-ff804cb656d7&ids=a72eec98-c5e0-4866-8953-765780acb6c1&ids=e77b2294-1bdd-4fba-928a-d81d2622312f&ids=965e01fc-318e-4c02-a801-d6fad60bfae4&ids=21ad5409-fe0b-4728-97e4-15520b9fc287&ids=1a777521-277c-4aeb-baf1-66871a7c2d2a&ids=c13a3449-9e0d-45a9-bcc0-518f55e45c8a&ids=5f2d329b-d59d-4112-b490-5114b830e34d&ids=bb966617-6c1f-4bb0-a1ed-ceb37ecade67&ids=05d11519-2b33-4742-aa87-3934632f2f2b&ids=39bfafe2-9628-434e-bd72-148051a47477&ids=481bea69-3cd5-45f3-8a52-2d4cc8fc8df7&ids=f95e407b-de69-416c-920c-6be8c9414862&ids=75940293-8fa6-47f9-ad5d-155b61933fdc&ids=e8e84ccf-f8a8-4551-9257-ef731d02116f&ids=e4991159-f088-4a2a-88b7-38d6ac47c6bc
```
```Shell
curl --remote-name --remote-header-name --request POST 'https://api.gdc.cancer.gov/data' --data @Payload
```
```Output
% Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
															 Dload  Upload   Total   Spent    Left  Speed
100 2563k    0 2562k  100   983   854k    327  0:00:03  0:00:03 --:--:--  776k
curl: Saved to filename 'gdc_download_20180830_132402.379282.tar.gz'
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

In this example we use `curl` to download a set of files from the GDC Legacy Archive; the payload is stored in a plain text file named `Payload`.


```Payload
{
    "ids":[
        "556e5e3f-0ab9-4b6c-aa62-c42f6a6cf20c",
        "e0de63e2-02f3-4309-9b24-69f4c24e85fc",
        "f1a06178-2ec2-4b06-83f3-3aedac332cfe",
        "11a8aca0-c8e6-4ff8-8ab6-fe18a1b8ba82",
        "69a69c84-00de-45ff-b397-fd2b6713ed4f",
        "9ec48233-395d-401e-b205-951c971f8dd4",
        "93129547-378c-4b69-b858-532abfff678e",
        "8d4277e9-a472-4590-886d-24dc2538ea65",
        "6733b412-56da-4f1c-a12b-ff804cb656d7",
        "a72eec98-c5e0-4866-8953-765780acb6c1",
        "e77b2294-1bdd-4fba-928a-d81d2622312f",
        "965e01fc-318e-4c02-a801-d6fad60bfae4",
        "21ad5409-fe0b-4728-97e4-15520b9fc287",
        "1a777521-277c-4aeb-baf1-66871a7c2d2a",
        "c13a3449-9e0d-45a9-bcc0-518f55e45c8a",
        "5f2d329b-d59d-4112-b490-5114b830e34d",
        "bb966617-6c1f-4bb0-a1ed-ceb37ecade67",
        "05d11519-2b33-4742-aa87-3934632f2f2b",
        "39bfafe2-9628-434e-bd72-148051a47477",
        "481bea69-3cd5-45f3-8a52-2d4cc8fc8df7",
        "f95e407b-de69-416c-920c-6be8c9414862",
        "75940293-8fa6-47f9-ad5d-155b61933fdc",
        "e8e84ccf-f8a8-4551-9257-ef731d02116f",
        "e4991159-f088-4a2a-88b7-38d6ac47c6bc"
    ]
}
```
```Shell
curl --remote-name --remote-header-name --request POST --header 'Content-Type: application/json' --data @request.txt 'https://api.gdc.cancer.gov/data'
```
```Output
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100 2562k    0 2561k  100  1145   788k    352  0:00:03  0:00:03 --:--:--  788k
curl: Saved to filename 'gdc_download_20160701_011007.tar.gz'
```

### Downloading Controlled-access Files

To download controlled-access files, a valid authentication token must be passed to the GDC API using the `X-Auth-Token` HTTP header:

```shell
token=$(<gdc-token-text-file.txt)

curl --remote-name --remote-header-name --header "X-Auth-Token: $token" 'https://api.gdc.cancer.gov/data/0eccf79d-1f1e-4205-910f-8e126b08276e'
```
```Output
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100 31.4M  100 31.4M    0     0   290k      0  0:01:50  0:01:50 --:--:--  172k
curl: Saved to filename 'ACOLD_p_TCGA_Batch17_SNP_N_GenomeWideSNP_6_A03_466078.tangent.copynumber.data.txt'
```

## Manifest endpoint

The `manifest` endpoint generates a download manifest file that can be used with the [GDC Data Transfer Tool](../../Data_Transfer_Tool/Users_Guide/Getting_Started.md). The Data Transfer Tool is recommended for transferring large volumes of data. The GDC API can also generate a download manifest from a list of results that match a [Search and Retrieval](Search_and_Retrieval.md) query. To do this, append `&return_type=manifest` to the end of the query.


### Using the manifest endpoint

The `manifest` endpoint allows users to create a download manifest, which can be used with the GDC Data Transfer Tool to download a large volume of data. The `manifest` endpoint generates a manifest file from a comma-separated list of UUIDs.

```shell
curl --remote-name --remote-header-name 'https://api.gdc.cancer.gov/v0/manifest/ae9db773-78ab-48d0-972d-debe1bedd37d,3d815e6e-db97-419d-ad7f-dba4e4023b3e'
```
```Output
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100   274  100   274    0     0   1042      0 --:--:-- --:--:-- --:--:--  1041
curl: Saved to filename 'gdc_manifest_20160428_234614.txt'
```

The `manifest` endpoint also supports HTTP POST requests in the same format as the `data` endpoint; see [above](#post-request-with-json-payload) for details.


### Using return_type=manifest

Alternatively, users can create a manifest by appending `&return_type=manifest` to a [Search and Retrieval](Search_and_Retrieval.md) query. In this example, we generate a download manifest for RNA-seq data files from solid tissue normal samples, that are part of the TCGA-KIRC project:

```Shell
curl --remote-name --remote-header-name 'https://api.gdc.cancer.gov/files?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B%22field%22%3A%22experimental_strategy%22%2C%22value%22%3A%5B%22RNA-Seq%22%5D%7D%7D%2C%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.project_id%22%2C%22value%22%3A%5B%22TCGA-KIRC%22%5D%7D%7D%2C%7B%22op%22%3A%22%3D%22%2C%22content%22%3A%7B%22field%22%3A%22cases.samples.sample_type%22%2C%22value%22%3A%5B%22Solid+Tissue+Normal%22%5D%7D%7D%5D%7D&size=30000&return_type=manifest'
```
```Output
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100 40663    0 40663    0     0  77109      0 --:--:-- --:--:-- --:--:-- 77306
curl: Saved to filename 'gdc_manifest.2016-06-28T13:26:33.850459.tsv'
```
