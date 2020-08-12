# Data Transfer Tool UI Release Notes


| Version | Date |
|---|---|
| [v0.6.0](DTT_UI_Release_Notes.md#v054) | August 12, 2020 |
| [v0.5.4](DTT_UI_Release_Notes.md#v054) | April 5, 2018 |
| [v0.5.3](DTT_UI_Release_Notes.md#v053) | December 14, 2017 |


## v0.5.4
* __GDC Product__: Data Transfer Tool UI
* __Release Date__: August 12, 2020


### New Features and Changes
* The server option can now be indicated within the interface. Previously the DTT-UI server defaulted to `api.gdc.cancer.gov`
* The DTT-UI now uses Data Transfer Tool v1.6, which uses Python3.  

### Bugs Fixed Since Last Release
* None

### Known Issues and Workarounds
* Download speeds for large numbers of small files may be better handled with the Command Line version of the Data Transfer Tool
* Data Submission to the GDC is not supported in the Data Transfer Tool UI.  Instead users must use the Command Line Data Transfer Tool

## v0.5.4
* __GDC Product__: Data Transfer Tool UI
* __Release Date__: April 5, 2018

### New Features and Changes
* None

### Bugs Fixed Since Last Release
* Download is now enabled for GDC reference and publication files. <!--SV-1045, DTT-100-->

### Known Issues and Workarounds
* Download speeds for large numbers of small files may be better handled with the Command Line version of the Data Transfer Tool
* Data Submission to the GDC is not supported in the Data Transfer Tool UI.  Instead users must use the Command Line Data Transfer Tool



## v0.5.3
* __GDC Product__: Data Transfer Tool UI
* __Release Date__: December 14, 2017


### New Features and Changes
* This is the first release for the Data Transfer Tool User Interface.  It allows users to download controlled access data using a simplified point and click interface.  This is a beta release and we welcome feedback on the user experience.  Important updates compared to the Command Line version include:
  * Upload and store authentication token between sessions
  * Easily view progress on a download manifest as the files are completed
  * View download history

### Bugs Fixed Since Last Release
* None

### Known Issues and Workarounds
* Download speeds for large numbers of small files may be better handled with the Command Line version of the Data Transfer Tool
* Data Submission to the GDC is not supported in the Data Transfer Tool UI.  Instead users must use the Command Line Data Transfer Tool
