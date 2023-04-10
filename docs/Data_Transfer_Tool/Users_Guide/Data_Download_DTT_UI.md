#Data Downloads with the Data Transfer Tool UI

##Data Transfer Tool UI: Overview
The UI version of the Data Transfer Tool was created for users who prefer a graphical interface over the command line or have limited command line experience.  The command line version is recommended for those users with more command line experience, require large data transfers of GDC data, or need to download a large numbers of data files.

###System Recommendations

The system recommendations for using the GDC Data Transfer Tool are as follows:

* OS: Linux (Ubuntu 14.x or later), OS X (10.9 Mavericks or later), or Windows (7 or later)
* CPU: At least four 64-bit cores, Intel or AMD
* RAM: At least 2 GiB
* Storage: Enterprise-class storage system capable of at least 1 Gb/s (gigabit per second) write throughput and sufficient free space for BAM files.

###Binary Distributions

Binary distributions are available on the GDC Transfer Tool page. To install the GDC Data Transfer UI download the respective binary distribution and unzip the distribution's archive to a location on the target system that is easily accessible.

###Binary Installation
Once the binary has been positioned in an appropriate location on the client's file system the application will need to run though a one-time installation process.  On first execution the binary install splash screen will appear showing the progress of the installation.  A hidden directory is created within the user's home directory labeled dtt that holds configuration and executable files.

![GDC DTT UI Installation](images/GDC_DTT_UI_INSTALLv7.png "GDC Data Transfer Tool UI Install")


###Preparing for Data Download

The GDC Data Transfer Tool UI is a stand-alone client application intended to work with data file information stored on the GDC Data Portals.  Data download information must first be gathered from the GDC Data Portal.  From there a manifest file can be [generated](https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Preparing_for_Data_Download_and_Upload/#obtaining-a-manifest-file-for-data-download) to supply the client.  Alternatively, individual file UUIDs can be provided to the UUID entry window located on the Download tab in the client.

![GDC DTT UI Start Page](images/DTT_UI_Start_Page.png)
##Downloads with UUIDs
The Data Transfer Tool UI can download files by individual UUID.  UUIDs can be entered into the client while on the download tab.  The single entry field labeled "Enter UUID(s)" allows the user to enter UUIDs individually.     

To obtain a data file's UUID from the GDC Data Portal, click on the file name to display the file's summary page which includes vital information such as its GDC UUID.  

##Downloads with Manifest
A portal-generated manifest file can be used with the Data Transfer Tool UI.  From the Download tab home page click on the Select Manifest File button.  A file system search window will popup allowing navigation to the manifest file.  

![GDC DTT UI Manifest Button Example](images/Manifest_button_DTT_UI_Start_Window.png "GDC Data Transfer Tool UI Manifest Button")     

##Download Progress Page  
The Download Progress Page is the command console for the Data Transfer Tool UI and allows users to monitor downloads. Progress of all downloads including the ability to start, stop, and restart a download are performed on the Download Progress Page.  Once file UUIDs or a manifest has been added to the queue the download can be started by clicking on the download button located at the lower right hand side of the page.  

![GDC DTT UI Download Progress Page_Download](images/Download_Progress_Page_download.png)

Once a download has completed, information about the downloads can be viewed from the Completed tab located at the bottom of the page.  Any Stopped or Failed downloads can also be viewed from their respective labels located at the bottom of the Status page.  

![DTT_UI_Download_Completed_Tab](images/DTT_UI_download_completed.png)

##Controlled Access File Downloads

Some files in the GDC are controlled access.  If you require access to these files please review the process outlined in the documentation [Obtaining Access to Controlled Data](https://gdc.cancer.gov/access-data/obtaining-access-controlled-data).    After appropriate authorization has been granted an access token can be generated to allow the Data Transfer Tool UI application access to the requested data files.  Documentation explaining the process of generating a token is located in the [Obtaining an Authentication Token for Data Downloads](https://docs.gdc.cancer.gov/Data_Transfer_Tool/Users_Guide/Preparing_for_Data_Download_and_Upload/#obtaining-an-authentication-token-for-data-downloads).  Once a token has been downloaded to a secure location on the client's local filesystem the Data Transfer Tool UI can now access it.  

![No_Token_icon](images/No_Token_file_dtt_ui.png)  

The current status of client authorization is viewable in the upper right corner of the application.  If the image and wording on the token manager access button is in red then no valid client token file has been uploaded.  To upload a valid token file click on the token status button.  The token manager window should appear allowing either a drag and drop token file upload or a file navigation window can be opened to navigate to the file location.

![Token Manager Window](images/Token_Manager_Window.png)

The token manager will verify access and display the projects for which the user has access.  To complete the token upload process click on the save button within the Token Manager window.

![Valid Token](images/validated_token.png)

##Settings and Advanced Settings

While the default download options will work for the majority of use cases, there are a vareity of ways to customize or modify the download process within the DTT UI.  Details of each of the settings are listed below.

![DTT Settings and Advanced Settings Page](images/DTT_Settings_Page.png)   

| Settings | Details
|----------|---------|
| Number of Client Connections: Default (3) | Number of concurrent client threads |
| Destination Folder: Default (User's Home Directory) | User selectable download file location |  
| Calculate Inbound Segment and check Md5sum on Restart: Default (On) | Verify previous partial downloaded files via segment check sum |
| Calculates check sums on previous downloaded files Default (On) | Verify downloaded files via file level check sum |
| Save Logs: - Download Navigation windows for client downloads | Export download or token log files via drop down and export log button |
| Debug Logging: Default (Off) | Enable debug level logging for file downloads |
| Block Size (Bytes): Default (1048576) | HTTP chunk size transfers |
| Save Interval (Bytes): Default (1000000) | save interval in bytes |
| Auto Retry: Default (On) | Enables auto retries of failed downloads |
| Retry(s): Default (5) | Number of retry attempts to download a file after failure |
| Seconds between Retrys: Default (5) | Number of seconds between retires |  
