#Data Downloads with the Data Transfer Tool UI

##Data Transfer Tool UI: Overview
The UI version of the Data Transfer Tool was created for users who prefer a graphical interface over the command line or have limited command line experience.  For users with more command line experience, require large data transfers of GDC data, or need to download a large numbers of data files the command line version is recommended.

###System Recommendations

The system recommendations for using the GDC Data Transfer Tool are as follows:

* OS: Linux (Ubuntu 14.x or later), OS X (10.9 Mavericks or later), or Windows (7 or later)
* CPU: At least four 64-bit cores, Intel or AMD
* RAM: At least 2 GiB
* Storage: Enterprise-class storage system capable of at least 1 Gb/s (gigabit per second) write throughput and sufficient free space for BAM files.

###Binary Distributions

Binary distributions are available on the GDC Transfer Tool page. To install the GDC Data Transfer UI Software download the respective binary distribution and unzip the distribution's archive to a location on the target system that is easily accessible.

###Binary Installation
Once the binary has been positioned in an appropriate location on the client's file system the application will need to run though a one time installation process.  On first execution the binary install splash screen will appear showing the progress of the installation.  A hidden directory is created within the user's home directory labeled dtt that holds configuration and executable files.

![GDC DTT UI Installation](images/GDC_DTT_UI_INSTALLv7.png "GDC Data Transfer Tool UI Install")


###Preparing for Data Download

The GDC Data Transfer Tool UI is a stand-alone client application intended to work with data file information stored on the GDC Data Portals.  Data download information must first be gathered from either the GDC Data Portal or Legacy Archive.  From there a manifest file can be [generated](https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Preparing_for_Data_Download_and_Upload/#obtaining-a-manifest-file-for-data-download) to supply the client.  Alternatively, individual file UUIDs can be provided to the UUID entry window located on the Download tab in the client.

![GDC DTT UI Start Page](images/DTT_UI_Start_Page.png)
##Downloads with UUIDs
The Data Transfer Tool UI can download files by individual UUID.  UUIDs can be entered into the client while on the download tab.  The single entry field labeled "Enter UUID(s)" allows the user to enter UUIDs individually.     

To obtain a data file's UUID from the GDC Data Portal, click on the file name to display the file's summary page which includes vital information such as its GDC UUID.  

##Downloads with Manifest
A portal-generated manifest file can be used with the Data Transfer Tool UI.  From the Download tab home page click on the Select Manifest File button.  A file system search window will popup allowing navigation to the manifest file.  

![GDC DTT UI Manifest Button Example](images/Manifest_button_DTT_UI_Start_Window.png "GDC Data Transfer Tool UI Manifest Button")     

##Download Progress Page  
The Data Transfer Tool Monitors downloads - The Download progress page is the command console for the Data Transfer Tool UI. Progress of all downloads including the ability to start and stop a download are performed on the Download Progress Page.   
