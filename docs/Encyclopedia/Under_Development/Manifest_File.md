# Manifest File #
## Description ##
A GDC manifest file is a file that is created by either the portal or API containing a list of files that the end user has requested.  
## Overview ##
When a user requests files from the GDC they are given the option to create a download manifest file that can be used in conjunction with the Data Transfer Tool to retrieve them.  The download manifest file contains the UUID, MD5 check sum, file size, and file name of each file listed in it.  The download manifest creation process is a function of the REST API and is made available to GDC users from an API download endpoint or from GDC portal's download cart.

For GDC submission projects that require multiple experimental data uploads a manifest file is required and needed for use with the Data Transfer Tool.  The format of the submission manifest differs from the download manifest because of the requirements needing to be met for upload.  A upload manifest file can be generated from either the API or submission portal.



## References ##
1. [Manifest endpoint](https://docs.gdc.cancer.gov/API/Users_Guide/Downloading_Files/#manifest-endpoint)
2. [GDC Data Transfer Tool](https://docs.gdc.cancer.gov/Data_Portal/Users_Guide/Cart/#gdc-data-transfer-tool)
3. [Upload Manifest](https://docs.gdc.cancer.gov/API/Users_Guide/Submission/#upload-manifest)
4. [Data Transfer Tool Submission with Manifest ](https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Preparing_for_Data_Download_and_Upload/#obtaining-a-manifest-file-for-data-uploads)
## External Links ##
* [external link name] (external link URL)

Categories: General
