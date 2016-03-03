# Data Transfer Tool Release Notes






## v0.3.20

* __GDC Product__: Data Transfer Tool
* __Release Date__: January 28, 2016


### New Features and Changes

* None to report

### Bugs Fixed Since Last Release

* GDC API compatibility fixes

### Known Issues and Workarounds

* On some terminals, dragging and dropping a file into the interactive client will add single quotes (' ') around the file path. This causes the interactive client to misinterpret the file path and generate an error when attempting to load a manifest file or token.
  * *Workaround:* Manually type out the file name or remove the single quotes from around the file path.

Release details are maintained in the [GDC Data Transfer Tool Change Log](https://github.com/NCI-GDC/gdc-client/blob/master/CHANGELOG.md)






## v0.2.18

* __GDC Product__: Data Transfer Tool
* __Release Date__: December 18, 2015

### New Features and Changes

* None to report

### Bugs Fixed Since Last Release

* Fixed a bug where adding a new line to the end of an auth token would break the header for a request.

### Known Issues and Workarounds

* On some terminals, dragging and dropping a file into the interactive client will add single quotes (' ') around the file path. This causes the interactive client to misinterpret the file path and generate an error when attempting to load a manifest file or token.
  * *Workaround:* Manually type out the file name or remove the single quotes from around the file path

Release details are maintained in the [GDC Data Transfer Tool Change Log](https://github.com/NCI-GDC/gdc-client/blob/master/CHANGELOG.md)






## v0.2.18-spr3

* __GDC Product__: Data Transfer Tool
* __Release Date__: November 18, 2015

### New Features and Changes

* Added molecular data upload via manifest files.
* Single part or multi-part upload is available dependent on file size.
* Users can resume partially uploaded files for multi-part uploads.

### Bugs Fixed Since Last Release

* None to report

### Known Issues and Workarounds

* On some terminals, dragging and dropping a file into the interactive client will add single quotes (' ') around the file path. This causes the interactive client to misinterpret the file path and generate an error when attempting to load a manifest file or token.
  * *Workaround:* Manually type out the file name or remove the single quotes from around the file path

Release details are maintained in the [GDC Data Transfer Tool Change Log](https://github.com/NCI-GDC/gdc-client/blob/master/CHANGELOG.md)






## v0.2.15.2

* __GDC Product__: Data Transfer Tool
* __Release Date__: August 31, 2015

### New Features and Changes

* None to report

### Bugs Fixed Since Last Release

* Fixed an issue where the authorization token was not properly passed to related files downloads, preventing download of controlled access related files.
* Fixed an issue where some error messages would not return a descriptive response for the user.
* Fixed an issue where the client would not properly read URLs if a trailing '/' was not added at the end.

### Known Issues and Workarounds

* On some terminals, dragging and dropping a file into the interactive client will add single quotes (' ') around the file path. This causes the interactive client to misinterpret the file path and generate an error when attempting to load a manifest file or token.
  * *Workaround:* Manually type out the file name or remove the single quotes from around the file path

Release details are maintained in the [GDC Data Transfer Tool Change Log](https://github.com/NCI-GDC/gdc-client/blob/master/CHANGELOG.md)






## v0.2.15.1

* __GDC Product__: Data Transfer Tool
* __Release Date__: August 7, 2015

### New Features and Changes

* Reorganized the High Performance Download Client into the GDC Client. This change was made in order to accomodate future functionality into one consolidated client.
* The GDC Client will now download associated files and annotations.

### Bugs Fixed Since Last Release

* None to report

### Known Issues and Workarounds

* On some Mac OS X systems, it is possible that a directory user permission error might be raised when attempting to download files.
  * *Workaround:* Switch the current working directory and attempt to download again

Release details are maintained in the [GDC Data Transfer Tool Change Log](https://github.com/NCI-GDC/gdc-client/blob/master/CHANGELOG.md)






## v0.2.13

* __GDC Product__: Data Transfer Tool
* __Release Date__: July 23, 2015

### New Features and Changes

* Removed debug definition from default configuration so users can specify prior to compile if they want verbose logging output

### Bugs Fixed Since Last Release

* Fixed a receive buffer space bug

### Known Issues and Workarounds

* None to report

Release details are maintained in the [GDC Data Transfer Tool Change Log](https://github.com/NCI-GDC/gdc-client/blob/master/CHANGELOG.md)






## v0.1.13

* __GDC Product__: Data Transfer Tool
* __Release Date__: June 2, 2015

### New Features and Changes

* Provides support for interactive text mode
* Adds a default server Uniform Resource Locator (URL)
* Updated command line syntax to default to TCP/HTTP in support of ease of use
* Provides support for Windows
* Packaged for OSX, Ubuntu, and Windows

### Bugs Fixed Since Last Release

* None to report

### Known Issues and Workarounds

* None to report

Release details are maintained in the [GDC Data Transfer Tool Change Log](https://github.com/NCI-GDC/gdc-client/blob/master/CHANGELOG.md)






## v0.1.10

* __GDC Product__: Data Transfer Tool
* __Release Date__: March 18, 2015

### New Features and Changes

* Allow users to retrieve large, high volume molecular data using a manifest file generated by the GDC Data Portal
* Allow users to download retrieved files

### Bugs Fixed Since Last Release

* None to report

### Known Issues and Workarounds

* Does not support the Windows platform

Release details are maintained in the [GDC Data Transfer Tool Change Log](https://github.com/NCI-GDC/gdc-client/blob/master/CHANGELOG.md)
