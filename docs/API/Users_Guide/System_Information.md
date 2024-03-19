# System Information

## Overview

The GDC API offers endpoints that provide information about the system. These endpoints are described below.


## GDC Notifications Endpoint

The `notifications` endpoint provides current user-facing notifications.

GDC notifications have a corresponding `level` with the following meanings:

| Level   | Meaning                                                              |
|---------|----------------------------------------------------------------------|
| INFO    | Non-essential information, e.g. regarding a new dataset              |
| WARNING | Important user information, e.g. regarding a dataset to be removed   |
| ERROR   | Important system information, e.g. regarding a GDC component         |
| DEBUG   | Unimportant system information, e.g. testing the notification system |

Notifications will indicate the GDC `components` to which they apply:

| Component   | Description                                                              |
|---------|----------------------------------------------------------------------|
| PORTAL   | The GDC Data Portal         |
| SUBMISSION   | The GDC Data Submission Portal |
| DOCUMENTATION | The GDC documentation site that contains GDC user guides, release notes, and the GDC Data Dictionary    |
| WEBSITE    | The GDC project website that includes information about the system. This does not include any of the above-listed GDC components.           |

### Sample Request

```Shell
curl --request GET https://api.gdc.cancer.gov/v0/notifications
```
```Response
{
  "data": [
    {
      "level": "INFO",
      "components": [
        "SUBMISSION_API"
        "LEGACY_API"
      ],
      "message": "The system is up!"
    }
  ]
}
```

## API Status Endpoint

The `status` endpoint provides information about the current status and version of the GDC API.

### Sample Request

``` shell
curl https://api.gdc.cancer.gov/status
```
``` python
import requests
import json

status_endpt = 'https://api.gdc.cancer.gov/status'
response = requests.get(status_endpt)
print json.dumps(response.json(), indent=2)
```
``` Reponse
{
  "commit": "74e1e3583c0f39fbf2149322addb7378206be3b9",
  "status": "OK",
  "tag": "1.2.0",
  "version": 1
}
```
