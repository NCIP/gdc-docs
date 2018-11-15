###Data Transfer Tool Configuration File  
The DTT has the ability to save configuration parameters in the format of a flat text file.


Example usage:

    gdc-client download d45ec02b-13c3-4afa-822d-443ccd3795ca --config my-dtt-config.dtt

Example of configuration file:

    [upload]
    path = /some/upload/path
    http_chunk_size = 1024


    [download]
    dir = /some/download/path
    http_chunk_size = 2048
    retry_amount = 6

Display Config Parameters:

    gdc-client settings download --config myconf.dtt
    [download]
    no_auto_retry = False
    no_file_md5sum = False
    save_interval = 1073741824
    http_chunk_size = 2048
    server = http://exmple-site.com
    n_processes = 8
    no_annotations = False
    no_related_files = False
    retry_amount = 6
    no_segment_md5sum = False
    manifest = []
    wait_time = 5.0
    no_verify = True
    dir = /some/download/path
