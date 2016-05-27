# Data Transfer Tool Release Notes





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

* Use of non-ASCII characters in token passed to Data Transfer Tool produces an incorrect error message.

### Known Issues and Workarounds

* On some terminals, dragging and dropping a file into the interactive client will add single quotes (' ') around the file path. This causes the interactive client to misinterpret the file path and generate an error when attempting to load a manifest file or token.
	* *Workaround:* Manually type out the file name or remove the single quotes from around the file path.
