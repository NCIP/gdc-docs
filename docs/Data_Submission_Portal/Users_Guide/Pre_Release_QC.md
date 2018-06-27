# Pre-Release Data Portal


## Getting Started


### The GDC Pre-Release Data Portal: An Overview

The Genomic Data Commons (GDC)  Portal provides users with web-based access to pre-released data from cancer genomics studies that have been harmonized by the GDC, but not yet released in the main GDC Data Portal. Key GDC Pre-Release Data Portal features include:

*   Access to data prior to release on the GDC Data Portal.
*   Repository page for browsing data by project / file / case
*   File / case faceted searches to filter data
*   Cart for collecting data files of interest
*   Authentication using eRA Commons credentials for access to controlled-access data files
*   Secure data download directly from the cart or using the [GDC Data Transfer Tool](https://gdc.cancer.gov/access-data/gdc-data-transfer-tool)
*   Use of API for query and download




## Navigation


Pre-Release Data Portal features are a subset of what can be found in the GDC Data Portal.  For more information on any of these general features please review the [Data Portal User Guide](/Data_Portal/Users_Guide/Getting_Started/#navigation).

[![GDC Views](images/AWG_Portal.png)](images/WG_Portal.png "Click to see the full image.")




## Authentication

### Overview

The GDC Pre-Release Data Portal provides access to datasets prior to release to a group of users specified by the data submitter.  This area is only available to data submitters (or their designees) for reviewing pre-release data.  Users must be granted access as specified in the admin portal section and also have downloader access within dbGaP for the specified project.

### GDC Authentication Tokens

The GDC Pre-Release Data Portal provides authentication tokens for use with the GDC Data Transfer Tool or the GDC API. To download a token:

1. Log into the GDC using your eRA Commons credentials
2. Click the username in the top right corner of the screen
3. Select the "Download token" option

![Token Download Button](images/gdc-data-portal-token-download.png)

A new token is generated each time the `Download Token` button is clicked.

For more information about authentication tokens, see [Data Security](../../Data/Data_Security/Data_Security.md#authentication-tokens).

**NOTE:** The authentication token should be kept in a secure location, as it allows access to all data accessible by the associated user account.

### Logging Out

To log out of the GDC, click the username in the top right corner of the screen, and select the Logout option.

![Logout link](images/gdc-data-portal-token-download.png)


## Data Transfer Tool

As with the GDC Data Portal, downloads of large or numerous files is best performed using the GDC Data Transfer Tool.  Information on the GDC Data Transfer Tool is available in the [GDC Data Transfer Tool User's Guide](/node/8196/).  An important distinction for use with the Pre-Release Data Portal is that it must always be used with a token and with the option `-s https://api.awg.gdc.cancer.gov`.

## GDC Pre-Release Data Admin Portal

### Overview

The GDC Pre-Release Data Admin Portal allows Pre-Release Data Portal admins to create and maintain Pre-Release Data Groups and associated projects, as well as grant appropriate access to users within these groups. To gain access to the Pre-Release Data Admin Portal please contact the GDC Helpdesk (support@nci-gdc.datacommons.io).

[![GDC Pre-Release Data Portal Main Page](images/AWG_Admin.png)](images/AWG_Admin.png "Click to see the full image.")

The Pre-Release Data Admin Portal is broken into two views on the left-most panel:

* __Users__: Allows admin to create, view, edit Pre-Release Data Portal user profiles
* __Groups__: Allows admin to manage groups projects / users

#### Definitions

| Entity | Definition |
|---|---|
| __User__  | An individual with an eRA Commons account. |
| __Project__  | A  collection of files and observations that are contained in the GDC database and have been registered in dbGAP as a project. Only certain projects are designated as Pre-Release Data projects.|
| __Group__  | A collection of users and projects.  When a user is assigned to a group, they will have access to the projects in that group when they login to the Pre-Release Data portal as long as they have downloader access to the project in dbGaP.|

### Users

The __Users__ section of the GDC Pre-Release Data Admin portal allows admins to manage and create Pre-Release Data users.

[![GDC Pre-Release Data Portal Main Page](images/AWG_Admin.png)](images/AWG_Admin.png "Click to see the full image.")

#### Creating Users

To create a new user in the Pre-Release Data Admin Portal, click on the `Create` button on the far right panel.

[![GDC Pre-Release Data Portal Main Page](images/AWG_Admin_Create_User.png)](images/AWG_Admin_Create_User.png "Click to see the full image.")

Then the following information must be supplied, before clicking the `Save` button:

* __eRA Commons ID__: The eRA Commons ID of the user to be added
* __Role__: Choose between `Admin` or `User` roles
* __Group (Optional)__: Choose existing groups to add the user to

After clicking `Save`, the user should appear in the list of users in the center panel.  Also clicking on the user in the list will display information about that user and gives the options to `Edit` the user profile, or `Delete` the user.

[![GDC Pre-Release Data Portal Main Page](images/AWG_Admin_New_User.png)](images/AWG_Admin_New_User.png "Click to see the full image.")

### Groups

The __Groups__ section of the GDC Pre-Release Data Admin portal allows admins to manage and create groups for which users and projects may be added.

[![GDC Pre-Release Data Portal Main Page](images/AWG_Admin_Group.png)](images/AWG_Admin_Group.png "Click to see the full image.")

#### Creating Groups

To create a new group in the Pre-Release Data Admin Portal, click on the `Create` button on the far right panel.

[![GDC Pre-Release Data Portal Main Page](images/AWG_Admin_Groups_Add.png)](images/AWG_Admin_Groups_Add.png "Click to see the full image.")

Then the following information must be supplied, before clicking the `Save` button:

* __Name__: The name of the group
* __Description__: The description of the group
* __Users (Optional)__: Choose existing users to add to the group
* __Projects(Optional)__: Choose existing projects to add to the group

After clicking `Save`, the group should appear in the list of groups in the center panel.  Also clicking on the group in the list will display information about that group and gives the options to `Edit` or `Delete` the group.

[![GDC Pre-Release Data Portal Main Page](images/AWG_Admin_New_Group.png)](images/AWG_Admin_New_Group.png "Click to see the full image.")

## API

API functionality is similar to what is available for the main GDC Data Portal.  You can read more about the GDC API in general in the [API User Guide](/API/Users_Guide/Getting_Started/).  Important differences for the AWG API include the following:

*  The base URL is different. Instead use https://api.awg.gdc.cancer.gov/
*  An authorization token must always be passed with every query rather than just for downloading controlled access data.
