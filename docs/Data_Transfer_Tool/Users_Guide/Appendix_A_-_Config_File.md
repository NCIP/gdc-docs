###Data Transfer Tool Configuration File  
The DTT has the ability to save and reuse configuration parameters in the format of a flat text file via a command line argument.  A simple text file needs to be created first with an extension of either txt or dtt.  The supported section headers are upload and download which can be used independently of each other or used in the same configuration file.  Each section header corresponds to the main functions of the application which are to either download data from the GDC portals or to upload data to the submission system of the GDC.  The configurable parameters are those listed in the help menus under either [download](https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Accessing_Built-in_Help/#download-help-menu) or [upload](https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Accessing_Built-in_Help/#upload-help-menu).            


Example usage:

    gdc-client download d45ec02b-13c3-4afa-822d-443ccd3795ca --config my-dtt-config.dtt

Example of configuration file:

    [upload]
    path = /some/upload/path
    upload_part_size = 1073741824


    [download]
    dir = /some/download/path
    http_chunk_size = 2048
    retry_amount = 6


###Display Config Parameters
This command line flag can be used with either the download or upload application feature to display what settings are active within a custom data transfer tool configuration file.  

    gdc-client settings download --config my-dtt-config.dtt
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
