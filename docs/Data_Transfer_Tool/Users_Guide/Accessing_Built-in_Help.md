

## Help Menus

The GDC Data Transfer Tool comes with built-in help menus. These menus are displayed when the GDC Data Transfer Tool is run with flags -h or --help for any of the main arguments to the tool. Running the GDC Data Transfer Tool without argument or flag will present a list of available command options.



```Shell
gdc-client --help
```
``` Output
usage: gdc-client [-h] [--version] {download,upload,interactive} ...

The Genomic Data Commons Command Line Client

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit

commands:
  {download,upload,interactive}
                        for more information, specify -h after a command
    download            download data from the GDC
    upload              upload data to the GDC
    interactive         run in interactive mode
```

 The available menus are provided below.

### Root menu

The GDC Data Transfer Tool displays the following help menu when executed without any arguments.

```Shell
gdc-client
```
```Output
usage: gdc-client [-h] [--version] {download,upload,interactive} ...
gdc-client: error: too few arguments
```


### Download help menu

The GDC Data Transfer Tool displays the following help menu for its download functionality.

```Shell
gdc-client download --help
```
```Output
usage: gdc-client download [-h] [--debug] [-v] [--log-file LOG_FILE]
                           [-T TOKEN | -t TOKEN] [-H HOST] [-P PORT] [-d DIR]
                           [-s server] [--no-segment-md5sums] [-n N_PROCESSES]
                           [--http-chunk-size HTTP_CHUNK_SIZE]
                           [--save-interval SAVE_INTERVAL]
                           [--no-related-files] [--no-annotations] [-u]
                           [-m MANIFEST]
                           [file_id [file_id ...]]

positional arguments:
  file_id               GDC files to download

optional arguments:
  -h, --help            show this help message and exit
  --debug               enable debug logging
  -v, --verbose         enable verbose logging
  --log-file LOG_FILE   log file [stderr]
  -t TOKEN, --token-file TOKEN
                        GDC API auth token file
  -H HOST, --host HOST  GDC API host [gdc-api.nci.nih.gov]
  -P PORT, --port PORT  GDC API port [443]
  -d DIR, --dir DIR     Directory to download files to. Defaults to current
                        dir
  -s server, --server server
                        The TCP server address server[:port]
  --no-segment-md5sums  Calculate inbound segment md5sums and/or verify
                        md5sums on restart
  -n N_PROCESSES, --n-processes N_PROCESSES
                        Number of client connections.
  --http-chunk-size HTTP_CHUNK_SIZE
                        Size in bytes of standard HTTP block size.
  --save-interval SAVE_INTERVAL
                        The number of chunks after which to flush state file.
                        A lower save interval will result in more frequent
                        printout but lower performance.
  --no-related-files    Do not download related files.
  --no-annotations      Do not download annotations.
  -u, --udt             Use the UDT protocol. Better for WAN connections
  -m MANIFEST, --manifest MANIFEST
                        GDC download manifest file
```

### Upload help menu

The GDC Data Transfer Tool displays the following help menu for its upload functionality.


```Shell
gdc-client upload --help
```
```Output
usage: gdc-client upload [-h] [--debug] [-v] [--log-file LOG_FILE]
                         [-T TOKEN | -t TOKEN] [-H HOST] [-P PORT]
                         [--project-id PROJECT_ID] [--identifier IDENTIFIER]
                         [--path path] [--upload-id UPLOAD_ID] [--insecure]
                         [--server SERVER] [--part-size PART_SIZE]
                         [-n N_PROCESSES] [--disable-multipart] [--abort]
                         [--resume] [--delete] [--manifest MANIFEST]

optional arguments:
  -h, --help            show this help message and exit
  --debug               enable debug logging
  -v, --verbose         enable verbose logging
  --log-file LOG_FILE   log file [stderr]
  -t TOKEN, --token-file TOKEN
                        GDC API auth token file
  -H HOST, --host HOST  GDC API host [gdc-api.nci.nih.gov]
  -P PORT, --port PORT  GDC API port [443]
  --project-id PROJECT_ID, -p PROJECT_ID
                        The project ID that owns the file
  --identifier IDENTIFIER, -i IDENTIFIER
                        The file id
  --path path, -f path  directory path to find file
  --upload-id UPLOAD_ID, -u UPLOAD_ID
                        Multipart upload id
  --insecure, -k        Allow connections to server without certs
  --server SERVER, -s SERVER
                        GDC API server address
  --part-size PART_SIZE, -ps PART_SIZE
                        Part size for multipart upload
  -n N_PROCESSES, --n-processes N_PROCESSES
                        Number of client connections
  --disable-multipart   Disable multipart upload
  --abort               Abort previous multipart upload
  --resume, -r          Resume previous multipart upload
  --delete              Delete an uploaded file
  --manifest MANIFEST, -m MANIFEST
                        Manifest which describes files to be uploaded
```
