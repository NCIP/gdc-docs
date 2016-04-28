# Downloading Files

GDC API's `data` endpoint allows users to download files stored in the GDC. Each file in GDC is assigned a Universally Unique Identifier (UUID). File UUIDs can be obtained from the GDC Data Portal in the form of individual UUIDs or as a manifest file that lists UUIDs of a group of files.

To download a file, users can pass individual UUID(s) to the `data` endpoint. If a single UUID is provided, the API will return the associated file. If a comma-separated list of UUIDs is provided, the API will return a single compressed (gzip) TAR file.

**Note:** Downloading controlled access data requires the use of a token. See [Authentication and Authorization](Authentication_and_Authorization.md).

## Downloading a Single File

To download a single file from the GDC, pass its UUID to the data endpoint.

```shell
curl -O -J 'https://gdc-api.nci.nih.gov/data/5b2974ad-f932-499b-90a3-93577a9f0573'
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
curl -O -J 'https://gdc-api.nci.nih.gov/data/e3228020-1c54-4521-9182-1ea14c5dc0f7,18e1e38e-0f0a-4a0e-918f-08e6201ea140'
```
```Output
% Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                               Dload  Upload   Total   Spent    Left  Speed
100  287k    0  287k    0     0  30131      0 --:--:--  0:00:09 --:--:-- 42759
curl: Saved to filename 'gdc_download_064d1aa8cc8cbab33e93979bebbf7d6af2d6a802.tar.gz'
```

## Downloading Controlled-access Files

To download controlled-access files, a valid authentication token must be passed to the GDC API using the `X-Auth-Token` HTTP header:

``` shell
export token=ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTOKEN-01234567890+AlPhAnUmErIcToKeN=0123456789-ALPHANUMERICTO

curl --remote-name --remote-header-name --header "X-Auth-Token: $token" 'https://gdc-api.nci.nih.gov/data/0eccf79d-1f1e-4205-910f-8e126b08276e'
```
```Output
  % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                 Dload  Upload   Total   Spent    Left  Speed
100 31.4M  100 31.4M    0     0   290k      0  0:01:50  0:01:50 --:--:--  172k
curl: Saved to filename 'ACOLD_p_TCGA_Batch17_SNP_N_GenomeWideSNP_6_A03_466078.tangent.copynumber.data.txt'
```
