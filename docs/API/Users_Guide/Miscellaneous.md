# Submission

## Overview

The GDC API also provides endpoints that are mainly used for programmatic access, but are documented here.

## GDC Notifications

The GDC API provides a `notifications` endpoint which returns any user facing notifications currently enabled.  GDC Notifications have a corresponding `level` with the following meanings.

| Level   | Meaning                                                              |
|---------+----------------------------------------------------------------------|
| INFO    | Non-essential information, e.g. regarding a new dataset              |
| WARNING | Important user information, e.g. regarding a dataset to be removed   |
| ERROR   | Important systen information, e.g. regarding a GDC component         |
| DEBUG   | Unimportant system information, e.g. testing the notification system |


```Command
curl --request GET https://gdc-api.nci.nih.gov/v0/notifications
```
```Response
{
  "data": [
    {
      "level": "INFO",
      "components": [
        "SUBMISSION_API",
        "LEGACY_API"
      ],
      "message": "The system is up!"
    }
  ]
}
```
