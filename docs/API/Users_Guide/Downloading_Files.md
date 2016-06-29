# Downloading Files

GDC API's `data` endpoint allows users to download files stored in the GDC by specifying file UUIDs.

To download a file, users can pass individual UUID(s) to the `data` endpoint. If a single UUID is provided, the API will return the associated file. If a comma-separated list of UUIDs is provided, the API will return an archive file containing the requested files.

The `manifest` endpoint generates a download manifest file that can be used with the GDC Data Transfer Tool. The Data Transfer Tool is recommended for transferring large volumes of data. The GDC API can also generate a download manifest from a list of results that match a [Search and Retrieval](Search_and_Retrieval.md) query. To do this, append `&return_type=manifest` to the end of the query.

**Note:** Downloading controlled access data requires the use of a token. See [Authentication and Authorization](Authentication_and_Authorization.md).

## Downloading a Single File

To download a single file from the GDC, pass its UUID to the data endpoint.

```shell
curl --remote-name --remote-header-name 'https://gdc-api.nci.nih.gov/data/5b2974ad-f932-499b-90a3-93577a9f0573'
```
```Output
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100 6111k  100 6111k    0     0   414k      0  0:00:14  0:00:14 --:--:--  412k
curl: Saved to filename '14-3-3_beta-R-V_GBL1112940.tif'
```

## Downloading Multiple Files

To download multiple files from the GDC, provide a comma-separated list of UUIDs to the `data` endpoint. The GDC API will return a `.tar.gz` archive containing the downloaded files.

```shell
curl --remote-name --remote-header-name 'https://gdc-api.nci.nih.gov/data/e3228020-1c54-4521-9182-1ea14c5dc0f7,18e1e38e-0f0a-4a0e-918f-08e6201ea140'
```
```Output
% Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                               Dload  Upload   Total   Spent    Left  Speed
100  287k    0  287k    0     0  30131      0 --:--:--  0:00:09 --:--:-- 42759
curl: Saved to filename 'gdc_download_064d1aa8cc8cbab33e93979bebbf7d6af2d6a802.tar.gz'
```

## Downloading Controlled-access Files

To download controlled-access files, a valid authentication token must be passed to the GDC API using the `X-Auth-Token` HTTP header:

```shell
export token=ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTO

curl --remote-name --remote-header-name --header "X-Auth-Token: $token" 'https://gdc-api.nci.nih.gov/data/0eccf79d-1f1e-4205-910f-8e126b08276e'
```
```Output
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100 31.4M  100 31.4M    0     0   290k      0  0:01:50  0:01:50 --:--:--  172k
curl: Saved to filename 'ACOLD_p_TCGA_Batch17_SNP_N_GenomeWideSNP_6_A03_466078.tangent.copynumber.data.txt'
```

## Creating a Manifest


### Using the manifest endpoint

The `manifest` endpoint allows users to create a download manifest, which can be used with the GDC Data Transfer Tool to download a large volume of data. The `manifest` endpoint generates a manifest file from a comma-separated list of UUIDs.

```shell
curl --remote-name --remote-header-name 'https://gdc-api.nci.nih.gov/v0/manifest/ae9db773-78ab-48d0-972d-debe1bedd37d,3d815e6e-db97-419d-ad7f-dba4e4023b3e'
```
```Output
% Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
															 Dload  Upload   Total   Spent    Left  Speed
100   274  100   274    0     0   1042      0 --:--:-- --:--:-- --:--:--  1041
curl: Saved to filename 'gdc_manifest_20160428_234614.txt'
```

### Using return_type=manifest

Alternatively, users can create a manifest by appending `&return_type=manifest` to a [Search and Retrieval](Search_and_Retrieval.md) query. In this example, we generate a download manifest for RNA-seq data files from solid tissue normal samples, that are part of the TCGA-KIRC project:

```Shell
curl --remote-name --remote-header-name 'https://gdc-api.nci.nih.gov/files?filters=%7B%0A%20%20%20%20%22op%22%3A%22and%22%2C%0A%20%20%20%20%22content%22%3A%5B%0A%20%20%20%20%20%20%20%20%7B%0A%20%20%20%20%20%20%20%20%20%20%20%20%22op%22%3A%22%3D%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%22content%22%3A%7B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22field%22%3A%22experimental_strategy%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22value%22%3A%5B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22RNA-Seq%22%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5D%0A%20%20%20%20%20%20%20%20%20%20%20%20%7D%0A%20%20%20%20%20%20%20%20%7D%2C%0A%20%20%20%20%20%20%20%20%7B%0A%20%20%20%20%20%20%20%20%20%20%20%20%22op%22%3A%22%3D%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%22content%22%3A%7B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22field%22%3A%22cases.project.project_id%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22value%22%3A%5B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22TCGA-KIRC%22%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5D%0A%20%20%20%20%20%20%20%20%20%20%20%20%7D%0A%20%20%20%20%20%20%20%20%7D%2C%0A%20%20%20%20%20%20%20%20%7B%0A%20%20%20%20%20%20%20%20%20%20%20%20%22op%22%3A%22%3D%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%22content%22%3A%7B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22field%22%3A%22cases.samples.sample_type%22%2C%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22value%22%3A%5B%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%22Solid%20Tissue%20Normal%22%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%5D%0A%20%20%20%20%20%20%20%20%20%20%20%20%7D%0A%20%20%20%20%20%20%20%20%7D%0A%20%20%20%20%5D%0A%7D&size=30000&return_type=manifest'
```
```Output
% Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
															 Dload  Upload   Total   Spent    Left  Speed
100 40663    0 40663    0     0  77109      0 --:--:-- --:--:-- --:--:-- 77306
curl: Saved to filename 'gdc_manifest.2016-06-28T13:26:33.850459.tsv'
```
