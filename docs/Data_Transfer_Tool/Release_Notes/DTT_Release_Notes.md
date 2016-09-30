# Data Transfer Tool Release Notes

## v1.1.0

* __GDC Product__: Data Transfer Tool
* __Release Date__: September 7, 2016


### New Features and Changes

* Partial extension added to all download files created during download. Removed after successful download.  
* Number of processes started by default changed to 8 (-n flag).

### Bugs Fixed Since Last Release

* None to report.

### Known Issues and Workarounds

* Use of non-ASCII characters in token passed to Data Transfer Tool will produce incorrect error message "Internal server error: Auth service temporarily unavailable".
* On some terminals, dragging and dropping a file into the interactive client will add single quotes (' ') around the file path. This causes the interactive client to misinterpret the file path and generate an error when attempting to load a manifest file or token.
	* *Workaround:* Manually type out the file name or remove the single quotes from around the file path.
<!---SV-457--->
* Use of a manifest file for uploads to the Submission Portal will produce an error message "ERROR: global name 'read_manifest' is not defined".
	* *Workaround:* Upload files via UUID instead or use the API/Submission Portal.



## v1.0.1

* __GDC Product__: Data Transfer Tool
* __Release Date__: June 2, 2016


### New Features and Changes

* MD5 checksum verification of downloaded files.
* BAM index files (.bai) are now automatically downloaded with parent BAM.
* UDT mode included to help improve certain high-speed transfers between the GDC and distant locations.

### Bugs Fixed Since Last Release

* None to report.

### Known Issues and Workarounds

* Use of non-ASCII characters in token passed to Data Transfer Tool will produce incorrect error message "Internal server error: Auth service temporarily unavailable".
* On some terminals, dragging and dropping a file into the interactive client will add single quotes (' ') around the file path. This causes the interactive client to misinterpret the file path and generate an error when attempting to load a manifest file or token.
	* *Workaround:* Manually type out the file name or remove the single quotes from around the file path.







## v1.0.0

* __GDC Product__: Data Transfer Tool
* __Release Date__: May 26, 2016


### New Features and Changes

* Single-thread and multi-threaded download capability
* User-friendly command line interface
* Progress bars provide visual representation of transfer status
* Optional interactive (REPL) mode
* Detailed help menus for upload and download functionality
* Support for authentication using a token file
* Support for authentication using a token string
* Resumption of incomplete uploads and downloads
* Initiation of transfers using manifests
* Initiation of transfers using file UUIDs
* Advanced configuration options
* Binary distributions available for Linux (Ubuntu), OS X, and Windows

### Bugs Fixed Since Last Release

* None to report.

### Known Issues and Workarounds

* Use of non-ASCII characters in token passed to Data Transfer Tool will produce incorrect error message "Internal server error: Auth service temporarily unavailable".
* On some terminals, dragging and dropping a file into the interactive client will add single quotes (' ') around the file path. This causes the interactive client to misinterpret the file path and generate an error when attempting to load a manifest file or token.
	* *Workaround:* Manually type out the file name or remove the single quotes from around the file path.
