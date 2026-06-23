

## Help Menus

The GDC Data Transfer Tool comes with built-in help menus. These menus are displayed when the GDC Data Transfer Tool is run with flags -h or --help for any of the main arguments to the tool. Running the GDC Data Transfer Tool without argument or flag will present a list of available command options.



```Shell
gdc-client --help
```
``` Output
usage: gdc-client [-h] [--version] {download,upload,settings} ...

The Genomic Data Commons Command Line Client

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit

commands:
  {download,upload,settings}
                        for more information, specify -h after a command
    download            download data from the GDC
    upload              upload data to the GDC
    settings            display default settings
```

 The available menus are provided below.

### Root menu

The GDC Data Transfer Tool displays the following output when executed without any arguments.

```Shell
gdc-client
```
```Output
usage: gdc-client [-h] [--version] {download,upload,settings} ...
gdc-client: error: too few arguments
```


### Download help menu

The GDC Data Transfer Tool displays the following help menu for its download functionality.

```Shell
gdc-client download --help
```
```Output
usage: gdc-client download [-h] [--debug]
                           [--log-file LOG_FILE]
                           [--color_off] [-t TOKEN_FILE]
                           [-d DIR] [-s server]
                           [--no-segment-md5sums]
                           [--no-file-md5sum]
                           [-n N_PROCESSES]
                           [--http-chunk-size HTTP_CHUNK_SIZE]
                           [--save-interval SAVE_INTERVAL]
                           [--no-verify]
                           [--no-related-files]
                           [--no-annotations]
                           [--no-auto-retry]
                           [--retry-amount RETRY_AMOUNT]
                           [--wait-time WAIT_TIME]
                           [--latest] [--config FILE] [-u]
                           [-m MANIFEST]
                           [file_id [file_id ...]]

positional arguments:
file_id               The GDC UUID of the file(s) to download

optional arguments:
-h, --help            show this help message and exit
--debug               Enable debug logging. If a failure occurs, the program
                      will stop.
--log-file LOG_FILE   Save logs to file. Amount logged affected by --debug
--color_off           Disable colored output
-t TOKEN_FILE, --token-file TOKEN_FILE
                      GDC API auth token file
-d DIR, --dir DIR     Directory to download files to. Defaults to current
                      dir
-s server, --server server
                      The TCP server address server[:port]
--no-segment-md5sums  Do not calculate inbound segment md5sums and/or do not
                      verify md5sums on restart
--no-file-md5sum      Do not verify file md5sum after download
-n N_PROCESSES, --n-processes N_PROCESSES
                      Number of client connections.
--http-chunk-size HTTP_CHUNK_SIZE, -c HTTP_CHUNK_SIZE
                      Size in bytes of standard HTTP block size.
--save-interval SAVE_INTERVAL
                      The number of chunks after which to flush state file.
                      A lower save interval will result in more frequent
                      printout but lower performance.
--no-verify           Perform insecure SSL connection and transfer
--no-related-files    Do not download related files.
--no-annotations      Do not download annotations.
--no-auto-retry       Ask before retrying to download a file
--retry-amount RETRY_AMOUNT
                      Number of times to retry a download
--wait-time WAIT_TIME
                      Amount of seconds to wait before retrying
--latest              Download latest version of a file if it exists
--config FILE         Path to INI-type config file
-u, --udt             Use the UDT protocol.
-m MANIFEST, --manifest MANIFEST
                      GDC download manifest file
```

### Upload help menu

The GDC Data Transfer Tool displays the following help menu for its upload functionality.


```Shell
gdc-client upload --help
```
```Output
usage: gdc-client upload [-h] [--debug]
                                            [--log-file LOG_FILE]
                                            [--color_off] [-t TOKEN_FILE]
                                            [--project-id PROJECT_ID]
                                            [--path path]
                                            [--upload-id UPLOAD_ID]
                                            [--insecure] [--server SERVER]
                                            [--part-size PART_SIZE]
                                            [--upload-part-size UPLOAD_PART_SIZE]
                                            [-n N_PROCESSES]
                                            [--disable-multipart] [--abort]
                                            [--resume] [--delete]
                                            [--manifest MANIFEST]
                                            [--config FILE]
                                            [file_id [file_id ...]]
positional arguments:
 file_id               The GDC UUID of the file(s) to upload

optional arguments:
 -h, --help            show this help message and exit
 --debug               Enable debug logging. If a failure occurs, the program
                       will stop.
 --log-file LOG_FILE   Save logs to file. Amount logged affected by --debug
 --color_off           Disable colored output
 -t TOKEN_FILE, --token-file TOKEN_FILE
                       GDC API auth token file
 --project-id PROJECT_ID, -p PROJECT_ID
                       The project ID that owns the file
 --path path, -f path  directory path to find file
 --upload-id UPLOAD_ID, -u UPLOAD_ID
                       Multipart upload id
 --insecure, -k        Allow connections to server without certs
 --server SERVER, -s SERVER
                       GDC API server address
 --part-size PART_SIZE
                       DEPRECATED in favor of [--upload-part-size]
 --upload-part-size UPLOAD_PART_SIZE, -c UPLOAD_PART_SIZE
                       Part size for multipart upload
 -n N_PROCESSES, --n-processes N_PROCESSES
                       Number of client connections
 --disable-multipart   Disable multipart upload
 --abort               Abort previous multipart upload
 --resume, -r          Resume previous multipart upload
 --delete              Delete an uploaded file
 --manifest MANIFEST, -m MANIFEST
                       Manifest which describes files to be uploaded
 --config FILE         Path to INI-type config file
```
