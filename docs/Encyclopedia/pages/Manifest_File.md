# Manifest File #
## Description ##
A GDC manifest file is a file that is created by either the portal or API containing a list of files that the end user has requested.  
## Overview ##
When a user requests files from the GDC they are given the option to create a download manifest file that can be used in conjunction with the Data Transfer Tool to retrieve them.  The download manifest file contains the UUID, MD5 checksum, file size, and file name of each file listed in it.  The download manifest creation process is a function of the GDC API and is made available to GDC users from an API download endpoint or from the GDC Data Portal.

A manifest file can also be used to simplify the upload process to the GDC Submission System.  The format of the submission manifest differs from the download manifest.  A upload manifest file can be generated from either the API or GDC Submission Portal.

## References ##
1. [Manifest endpoint](/API/Users_Guide/Downloading_Files/#manifest-endpoint)
2. [GDC Data Transfer Tool](/Data_Portal/Users_Guide/Cart/#gdc-data-transfer-tool)
3. [Upload Manifest](/API/Users_Guide/Submission/#upload-manifest)
4. [Data Transfer Tool Submission with Manifest ](/Data_Transfer_Tool/Users_Guide/Preparing_for_Data_Download_and_Upload/#obtaining-a-manifest-file-for-data-uploads)

## External Links ##
* N/A

Categories: General
