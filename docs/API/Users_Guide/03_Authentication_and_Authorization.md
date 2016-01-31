## 3.1 Overview

The GDC API and allows queries on available data and download of open access data without authentication. To download controlled access data, users are required to have a NIH eRA Commons account and appropriate dbGaP access. As described in the next section, these credentials can be used to obtain a API token via the GDC Portal.

## 3.2 Obtaining a Token
<a name="obtaining-a-token"></a>

Tokens are strings of characters provided with every API query for authentication. A token can be obtained through the GDC Data Portal (or GDC Data Submission Portal) once the user's NIH eRA Commons account has been established. Your GDC account will be automatically established upon first login and the user will see a drop down menu next to her/his username.

[[nid:8265]]

Image 3.1-1 - GDC Login and Download Token Dropdown

Users can download and save your token from here.

## 3.3 Using the Token

Each GDC API request must include a "X-Auth-Token" custom header.

**Example**:

    export token=YOUR_TOKEN
    curl -H 'X-Auth-Token:$token'  "https://gdcapi.nci.nih.gov/data/49ac8944-1468-456a-bb65-b08c7e24a97a"

In the example above, replace YOUR_TOKEN with the token downloaded from the portal.

## 3.4 Token Expiration

Tokens are valid for ninety days from the time of download. Using an expired token will result in a 401 HTTP error code:

    HTTP/1.1 403 FORBIDDEN{
      "error": "You don't have access to the data"
    }

**NOTE**: Invalid credentials will result in a server error even if the resource is open access.
