# Downloading Files

GDC API's `data` endpoint allows users to download files stored in the GDC. Each file in GDC is assigned a Universally Unique Identifier (UUID). File UUIDs can be obtained from the GDC Data Portal in the form of individual UUIDs or as a manifest file that lists UUIDs of a group of files.

To download a file, users can pass individual UUID(s) to the `data` endpoint. If a single UUID is provided, the API will return the associated file. If a comma-separated list of UUIDs is provided, the API will return a single compressed (gzip) TAR file.

**Note:** Downloading controlled access data requires the use of a token. See [Authentication and Authorization](Authentication_and_Authorization.md).

## Downloading a Single File

To download a single file from the GDC, pass its UUID to the data endpoint.

```shell
curl -O -J 'https://gdc-api.nci.nih.gov/data/96487cd7-8fa8-4bee-9863-17004a70b2e9'
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
curl -O -J 'https://gdc-api.nci.nih.gov/data/96487cd7-8fa8-4bee-9863-17004a70b2e9,5e55748f-61fa-43e8-886f-6f1ec7a91af6'
```
```Output
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100 10.1M    0 10.1M    0     0   341k      0 --:--:--  0:00:30 --:--:--  283k
curl: Saved to filename 'gdc_download_cfbc1fe89d423c5dc60b5eecf12b797d14553d17.tar.gz'
```

## Downloading Controlled-access Files

To download controlled-access files, a valid authentication token must be passed to the GDC API using the `X-Auth-Token` HTTP header:

``` shell
export token=ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTO

curl -O -J -H "X-Auth-Token: $token" 'https://gdc-api.nci.nih.gov/data/a1c1b23b-cc41-4e85-b1b7-62a42873c5af'
```
```Output
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100 31.4M  100 31.4M    0     0   290k      0  0:01:50  0:01:50 --:--:--  172k
curl: Saved to filename 'ACOLD_p_TCGA_Batch17_SNP_N_GenomeWideSNP_6_A03_466078.tangent.copynumber.data.txt'
```
