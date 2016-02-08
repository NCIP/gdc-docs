# 3 - Accessing Built-in Help

## Help Menus

The GDC Data Transfer Tool comes with built-in help menus. These menus are displayed when the GDC Data Transfer Tool is run without commands or when -h / --help flag is specified for any command. Running the GDC Data Transfer Tool without commands will present a list of available commands. The available menus are provided below.

### The GDC Data Transfer Tool Built-in Help Menus (run without commands)

    > gdc-client
    Type 'help' for a list of commands or 'help <topic>' for detailed usage.
    
    Basic commands are:
    - download   (download files in registry)
    - add        (adds ids to registry)
    - list       (lists file ids already registered)
    - manifest   (add ids from a GDC manifest file to registry)
    - remove     (remove ids from registry)
    - token      (load an authorization token file)
    - cd         (move to directory you want to download to)
    - pwd        (print the current working directory)
    - set        (set advanced configuration setting)
    - settings   (list advanced configuration settings)
    - upload     (upload files to object storage)
    - delete     (delete files from object storage)
    - abort      (abort a previous partial upload)
    
    TIPS:
    - Rather than type out path names, try dragging and dropping manifest and token files into the terminal.
    - You can execute shell commands by prepending '!', i.e. !ls.
    - You can run the gdc-client binary with advanced options from the command line (gdc-client --help).
    
    gdc-client repl >


### GDC Data Transfer Tool Download Help Menus (with -h / --help flag specified)

    > gdc-client download -h
    usage: gdc-client download [-h] [-m MANIFEST] [-v] [-d DIR] [-s server]
                               [--no-segment-md5sums] [--debug] [-n N_PROCESSES]
                               [--http-chunk-size HTTP_CHUNK_SIZE]
                               [--save-interval SAVE_INTERVAL]
                               [--no-related-files] [--no-annotations]
                               [-t TOKEN | -T TOKEN] [-u] [-H PROXY_HOST]
                               [-P PROXY_PORT] [-e]
                               [file_id [file_id ...]]
    
    positional arguments:
      file_id               uuids to download
    
    optional arguments:
      -h, --help            show this help message and exit
      -m MANIFEST, --manifest MANIFEST
                            GDC Download manifest file.
      -v, --verbose         verbose logging
      -d DIR, --dir DIR     Directory to download files to. Defaults to current
                            dir
      -s server, --server server
                            The UDT server address server[:port]
      --no-segment-md5sums  Calculate inbound segment md5sums and/or verify
                            md5sums on restart
      --debug               Print stack traces
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
      -t TOKEN, --token-file TOKEN
                            authentication token file
      -T TOKEN, --token TOKEN
                            authentication token
      -u, --udt             Use the UDT protocol. Better for WAN connections
      -H PROXY_HOST, --proxy-host PROXY_HOST
                            The port to bind the local proxy to
      -P PROXY_PORT, --proxy-port PROXY_PORT
                            The port to bind the local proxy to
      -e, --external-proxy  Do not create a local proxy but bind to an external
                            one
    >


### GDC Data Transfer Tool Upload Help Menus (with -h / --help flag specified)

    > gdc-client upload -h
    usage: gdc-client upload [-h] [--project-id PROJECT_ID]
                             [--identifier IDENTIFIER] [--path path] --token file
                             [--insecure] [--verbose] [--server SERVER]
                             [--part-size PART_SIZE] [-n N_PROCESSES]
                             [--upload-id UPLOAD_ID] [--disable-multipart]
                             [--abort] [--resume] [--delete] [--manifest MANIFEST]
    
    optional arguments:
      -h, --help            show this help message and exit
      --project-id PROJECT_ID, -p PROJECT_ID
                            The project ID that owns the file
      --identifier IDENTIFIER, -i IDENTIFIER
                            The id or alias
      --path path, -f path  directory path to find file
      --token file, -t file
                            auth token
      --insecure, -k        Allow connections to server without certs
      --verbose, -v         Print stack traces
      --server SERVER, -s SERVER
                            GDC API server address
      --part-size PART_SIZE, -ps PART_SIZE
                            Part size for multipart upload
      -n N_PROCESSES, --n-processes N_PROCESSES
                            Number of client connections.
      --upload-id UPLOAD_ID, -u UPLOAD_ID
                            Multipart upload id
      --disable-multipart   Disable multipart upload
      --abort               Abort previous multipart upload
      --resume, -r          Resume previous multipart upload
      --delete              Delete an uploaded file
      --manifest MANIFEST, -m MANIFEST
                            Manifest which describes files to be uploaded

