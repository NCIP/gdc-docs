# Controlled Data Access
The GDC API and allows queries on available data and download of open access data without authentication. To download controlled access data, users are required to have a NIH eRA Commons account and appropriate dbGaP access. Clink on these links to learn more about how to <a href="https://gdc.nci.nih.gov/access-data/obtaining-access-controlled-data" target="_blank">obtain access to controlled data</a> or how to <a href="https://gdc.nci.nih.gov/submit-data/obtaining-access-submit-data" target="_blank">submit data</a> to GDC.

As described in the following section, these credentials are used to obtain a GDC API token via the GDC Portals.

## Obtaining a Token
Tokens are strings of characters provided with every API query for authentication. A token can be obtained through the <a href="https://gdc-portal.nci.nih.gov" target="_blank">GDC Data Portal</a> or <a href="https://gdc-portal.nci.nih.gov/submission" target="_blank">GDC Data Submission Portal</a> once the user's NIH eRA Commons account has been established. Your GDC account will be automatically established upon first login and the user will see a drop down menu next to her/his username.

![GDC Login and Download Token Dropdown Menu](images/03-01__GDC_Login_and_Download_Token_Dropdown_Menu.png)  
**GDC Login and Download Token Dropdown Menu**

Users can download and save their token from here.

## Using the Token
Each GDC API request must include a "X-Auth-Token" custom header.

```
    export token=YOUR_TOKEN
    curl -H 'X-Auth-Token:$token'  "https://gdcapi.nci.nih.gov/data/49ac8944-1468-456a-bb65-b08c7e24a97a"
```

>In the example above, replace **YOUR_TOKEN** with the token downloaded from the portal.

## Token Expiration
Tokens are valid for 90 days from the time of download. Please note:

- Using an expired token will result in a 401 HTTP error code.
- Even if the resource is open access, using invalid credentials will result in a server error.

```
    HTTP/1.1 403 FORBIDDEN{
      "error": "You don't have access to the data"
    }
```

**NOTE**: 
