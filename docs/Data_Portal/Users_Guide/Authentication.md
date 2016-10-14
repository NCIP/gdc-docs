# Authentication

## Overview

The GDC Data Portal provides granular metadata for all datasets available in the GDC. Any user can see a listing of all available data files, including controlled-access files. The GDC Data Portal also allows users to download open-access files without logging in, but downloading of controlled-access files is restricted to authorized users and requires authentication.

Controlled-access files are identified using a "lock" icon:

[![GDC Data Portal Main Page](images/gdc-data-portal-controlled-files.png)](images/gdc-data-portal-controlled-files.png "Click to see the full image.")

Upon successful authentication, GDC Data Portal users can:

- see which controlled-access files they have access to;
- download controlled-access files directly from the GDC Data Portal;
- download an authentication token for use with the GDC Data Transfer Tool or the GDC API.

The rest of this section describes controlled data access features of the GDC Data Portal available to authorized users. For more information about open and controlled-access data, and about obtaining access to controlled data, see [Data Access Processes and Tools](https://gdc.cancer.gov/access-data/data-access-processes-and-tools).

## Authentication via eRA Commons

Users can log into the GDC using their eRA Commons credentials by clicking the "Login" button in the upper right corner of the screen. If authentication is successful, the eRA Commons username will be displayed in the upper right corner of the screen, in place of the "Login" button.

## GDC Authentication Tokens

The GDC Data Portal provides authentication tokens for use with the GDC Data Transfer Tool or the GDC API. To download a token:

0. Log into the GDC using your eRA Commons credentials,
0. Click the username in the top right corner of the screen,
0. Select the "Download token" option.

![Token Download Button](images/gdc-data-portal-token-download.png)

For more information about authentication tokens, see [Data Security](../../Data/Data_Security/Data_Security.md#authentication-tokens).

**NOTE:** The authentication token should be kept in a secure location, as it allows access to all data accessible by the associated user account.

## Display Only Authorized Datasets

The "Only My Projects" checkbox, available in [Data view](Data_View.md) to authenticated users, limits data file listings to only those files to which the user has access.

[![Only My Projects checkbox](images/gdc-data-portal-only-my-projects.png)](images/gdc-data-portal-only-my-projects.png "Click to see the full image.")


## Viewing File Authorization Information

The GDC Data Portal displays the following file authorization information for authenticated users:

- Authorization status of individual files in the "My Projects" column, available in [Data view](Data_View.md) and [cart](Cart.md) file listings.
- Summary of file authorization status in the user's cart.

[![File Authorization Summary in GDC Cart](images/gdc-portal-cart-authorization-summary.png)](images/gdc-portal-cart-authorization-summary.png "Click to see the full image.")



## Logging Out

To log out of the GDC, click the username in the top right corner of the screen, and select the Logout option.

![Logout link](images/gdc-data-portal-token-download.png)
