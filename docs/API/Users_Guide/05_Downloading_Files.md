## 5.1 Overview

This section explains how to use a GDC API endpoint to download single or multiple files. Each file in GDC is assigned a Universally Unique Identifier (UUID). The UUIDs of single or multiple files are needed to download them. UUIDs can be obtained by downloading a manifest from the Data Portal, or directly from the file entity page on the Data Portal. If a single UUID is provided, the API will return just the associated file. If a list of UUIDs are provided, the API will package all of the associated files together and return a single compressed (gzip) TAR file.

The GDC Data Portal uses the GDC API endpoint allows users download files via a web browser. However, the download endpoint can be called directly from command tools like curl. For example, the basic syntax for performing a download using the curl command is below:

    curl -J -O -L https://gdc-api.nci.nih.gov/data/UUID

**Note:** Downloading controlled access data requires the use of a token. See section 3 - Authentication and Authorization.

## 5.2 Downloading Files

### 5.2.1 Downloading a Single File

To download a single file from GDC, simply send its UUID to the data endpoint.

    curl --O -J -L 'https://gdc-api.nci.nih.gov/data/acad3917-f42e-40f9-b9dc-f6cca2909273'
      % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                     Dload  Upload   Total   Spent    Left  Speed
    100 7905k  100 7905k    0     0  7549k      0  0:00:01  0:00:01 --:--:-- 7550k
    curl: Saved to filename '3999492009_R01C01_Grn.idat'

#### 5.2.1.1 Avoid Duplicate Names

Although all files stored within the GDC are assigned a GDC UUID, the actual name of the file may not always be unique. For example, there are many files named 'README.txt' associated with different file archives imported from the The Cancer Genome Atlas (TCGA) program.

The following example using various curl ensures the downloaded file is stored locally using a unique file name.
```
    $ #  assign file uuid to a variable
    $ uuid=acad3917-f42e-40f9-b9dc-f6cca2909273
    $
    $ # url variable
    $ url=https://gdc-api.nci.nih.gov/data/$uuid
    $
    $ echo $url
    https://gdc-api.nci.nih.gov/data/acad3917-f42e-40f9-b9dc-f6cca2909273
    $
    $ # get remote filename by parsing  http header returned by curl -I option
    $ filename=$(curl -I  $url | grep -o -E 'filename=.*$' | tr -d '\r\n' | sed -e 's/filename=//')
      % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                     Dload  Upload   Total   Spent    Left  Speed
      0 7905k    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0
    $ echo $filename
    3999492009_R01C01_Grn.idat
    $
    $ #  append the uuid to remote filename to create unique filename
    $ filename=$uuid.$filename
    $ echo $filename
    acad3917-f42e-40f9-b9dc-f6cca2909273.3999492009_R01C01_Grn.idat
    $
    $ # rename the downloaded filename using the curl -o option
    $ curl -o $filename -L $url
      % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                     Dload  Upload   Total   Spent    Left  Speed
    100 7905k  100 7905k    0     0  6608k      0  0:00:01  0:00:01 --:--:-- 6615k
    $
    $ ls -l $filename
    -rw-rw-r-- 1 ubuntu ubuntu 8095289 Mar 19 20:39 acad3917-f42e-40f9-b9dc-f6cca2909273.3999492009_R01C01_Grn.idat
    $
```

### 5.2.2 Downloading Multiple Files

If a comma-separated list of UUIDs is provided to the data endpoint, the GDC REST API returns a uniquely named tar.gz compressed archive file containing the downloaded files.

```
    $ curl -O -J -L  'https://gdc-api.nci.nih.gov/data/acad3917-f42e-40f9-b9dc-f6cca2909273,8cfdebb1-91d4-48e4-a6b0-e2c71041a13b'
      % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
                                     Dload  Upload   Total   Spent    Left  Speed
    100 8123k    0 8123k    0     0  2827k      0 --:--:--  0:00:02 --:--:-- 2827k
    curl: Saved to filename 'gdc_download_5a40a579009104cabec3b1cf656ce57e8f2d8ec6.tar.gz'
```
