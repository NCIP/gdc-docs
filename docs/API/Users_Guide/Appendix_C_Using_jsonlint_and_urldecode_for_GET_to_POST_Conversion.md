# Appendix C - Using jsonlint and urldecode for GET to POST Conversion

## Example: GET to POST conversion using jsonlint and urldecode

The tools described in section 2.2 of this user's guide convert an API query from GET to POST and ensure the query is correct before submitting it to the API.

The first operation consists in breaking up the query below by identifying its parameters

    https://gdc-api.nci.nih.gov/files?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.clinical.gender%22%2C%22value%22%3A%5B%22female%22%5D%7D%7D%2C%7B%22op%22%3A%22%3E%3D%22%2C%22content%22%3A%7B%22field%22%3A%22cases.clinical.age_at_diagnosis%22%2C%22value%22%3A%5B%229125%22%5D%7D%7D%2C%7B%22op%22%3A%22%3C%3D%22%2C%22content%22%3A%7B%22field%22%3A%22cases.clinical.age_at_diagnosis%22%2C%22value%22%3A%5B%2212775%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.primary_site%22%2C%22value%22%3A%5B%22Skin%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_format%22%2C%22value%22%3A%5B%22CEL%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.case_id%22%2C%22value%22%3A%5B%22ccf8809a-cdfd-46f7-b607-001b8b91c1a4%22%5D%7D%7D%5D%7D&from=1&size=2&fields=file_id,file_size&pretty=true

The following parameters can be identified: 

* filters 
* from 
* size 
* fields

**Note**: the pretty paremeter can be discarded at this point, since tools used for POST usually provide JSON parsing.

The resulting skeleton for the POST query would be:

    {
      "from": 1,
      "size": 2,
      "fields": "file_id,file_size",
      "filters": "TBD"
    }

**Note**: Filters being a complex string it has been replaced in the above skeleton by "TBD", this is sufficient to validate the skeleton.

At this point, [JSONLint](http://jsonlint.com/) can be used to ensure the skeleton does not contain errors. Copying the skeleton in [JSONLint](http://jsonlint.com/) query string and clicking "Validate" should return "Valid JSON" at the bottom of the page.

The "filters" parameter contains a complex string:

    %7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.clinical.gender%22%2C%22value%22%3A%5B%22female%22%5D%7D%7D%2C%7B%22op%22%3A%22%3E%3D%22%2C%22content%22%3A%7B%22field%22%3A%22cases.clinical.age_at_diagnosis%22%2C%22value%22%3A%5B%229125%22%5D%7D%7D%2C%7B%22op%22%3A%22%3C%3D%22%2C%22content%22%3A%7B%22field%22%3A%22cases.clinical.age_at_diagnosis%22%2C%22value%22%3A%5B%2212775%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.primary_site%22%2C%22value%22%3A%5B%22Skin%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_format%22%2C%22value%22%3A%5B%22CEL%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.case_id%22%2C%22value%22%3A%5B%22ccf8809a-cdfd-46f7-b607-001b8b91c1a4%22%5D%7D%7D%5D%7D

[UrlDecode.org](http://urldecode.org/) provides the ability to decode this string into a JSON.

    {"op":"and","content":[{"op":"in","content":{"field":"cases.clinical.gender","value":["female"]}},{"op":">=","content":{"field":"cases.clinical.age_at_diagnosis","value":["9125"]}},{"op":"<=","content":{"field":"cases.clinical.age_at_diagnosis","value":["12775"]}},{"op":"in","content":{"field":"cases.project.primary_site","value":["Skin"]}},{"op":"in","content":{"field":"files.data_format","value":["CEL"]}},{"op":"in","content":{"field":"cases.case_id","value":["ccf8809a-cdfd-46f7-b607-001b8b91c1a4"]}}]}

The last action consists in replacing "TBD" in the skeleton with the decoded JSON.

    {
      "from": 1,
      "size": 2,
      "fields": "file_id,file_size",
      "filters": {"op":"and","content":[{"op":"in","content":{"field":"cases.clinical.gender","value":["female"]}},{"op":">=","content":{"field":"cases.clinical.age_at_diagnosis","value":["9125"]}},{"op":"<=","content":{"field":"cases.clinical.age_at_diagnosis","value":["12775"]}},{"op":"in","content":{"field":"cases.project.primary_site","value":["Skin"]}},{"op":"in","content":{"field":"files.data_format","value":["CEL"]}},{"op":"in","content":{"field":"cases.case_id","value":["ccf8809a-cdfd-46f7-b607-001b8b91c1a4"]}}]}
    }

If desired, [JSONLint](http://jsonlint.com/) can be used to prettify the above JSON.

    //https://gdc-api.nci.nih.gov/files
    {
      "from": 1,
      "size": 2,
      "fields": "file_id,file_size",
      "filters": {
        "op": "and",
        "content": [
          {
            "op": "in",
            "content": {
              "field": "cases.clinical.gender",
              "value": [
                "female"
              ]
            }
          },
          {
            "op": ">=",
            "content": {
              "field": "cases.clinical.age_at_diagnosis",
              "value": [
                "9125"
              ]
            }
          },
          {
            "op": "<=",
            "content": {
              "field": "cases.clinical.age_at_diagnosis",
              "value": [
                "12775"
              ]
            }
          },
          {
            "op": "in",
            "content": {
              "field": "cases.project.primary_site",
              "value": [
                "Skin"
              ]
            }
          },
          {
            "op": "in",
            "content": {
              "field": "files.data_format",
              "value": [
                "CEL"
              ]
            }
          },
          {
            "op": "in",
            "content": {
              "field": "cases.case_id",
              "value": [
                "ccf8809a-cdfd-46f7-b607-001b8b91c1a4"
              ]
            }
          }
        ]
      }
    }
