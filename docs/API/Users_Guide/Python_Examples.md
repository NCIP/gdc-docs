# Using Python to Access the GDC API

Python can be a versatile tool for retrieving information from the GDC API and performing downstream processing. This guide details some examples that demonstrate the basic methods for querying the API using Python. Most of the examples in this guide will use the [requests](http://docs.python-requests.org/en/master/) Python library and should be compatible with Python3.

## Querying metadata

### A Basic Query

This example passes some basic parameters (fields, format, size) to the `cases` endpoint and prints the results.

```Python
import requests
import json

cases_endpt = 'https://api.gdc.cancer.gov/cases'

fields = ["submitter_id","case_id","primary_site","disease_type","diagnoses.vital_status"]
fields = ','.join(fields)

params = {
    "fields":fields,
    "format":"TSV",
    "size":"100",
    }

response = requests.get(cases_endpt, params = params)

print(response.content)
```

### A Filtered Query

In the next example, a `filters` parameter is added to the query. This parameter is passed as a Python dictionary object. The filter used in this example will only show cases that come from a patient with kidney disease (`primary_site: Kidney`).

```Python
import requests
import json


fields = ["submitter_id","case_id","primary_site","disease_type","diagnoses.vital_status"]
fields = ','.join(fields)

cases_endpt = 'https://api.gdc.cancer.gov/cases'

filters = {
    "op":"in",
    "content":{
        "field":"primary_site",
        "value":["Kidney"]
            }
        }

# The filters parameter needs to be converted from a dictionary to JSON object

params = {
    "filters":json.dumps(filters),
    "fields":fields,
    "format":"TSV",
    "size":"100",
    }

response = requests.get(cases_endpt, params = params)

print(response.content)
```

### A More Complicated Filter

The following example returns information about the RNA-Seq BAM files that originate from lung cancer patients.   


```Python
import requests
import json

fields = ["file_name","cases.submitter_id","cases.samples.sample_type",
    "cases.disease_type","cases.project.project_id"]

fields = ','.join(fields)

files_endpt = 'https://api.gdc.cancer.gov/files'

filters = {
    "op":"and",
    "content":[
        {
        "op":"in",
        "content":{
            "field":"cases.project.primary_site",
            "value":["Lung"]
            }
        },
        {
        "op":"in",
        "content":{
            "field":"files.experimental_strategy",
            "value":["RNA-Seq"]
            }
        },
        {
        "op":"in",
        "content":{
            "field":"files.data_format",
            "value":["BAM"]
            }
        }
    ]
}

params = {
    "filters":json.dumps(filters),
    "fields":fields,
    "format":"TSV",
    "size":"2000"
    }

response = requests.get(files_endpt, params = params)

print(response.content)
```

## Troubleshooting

The following script should produce an unformatted JSON string when run. Run this to verify that a valid connection is being made to the GDC API.  

```Python
import requests
files_endpt = 'https://api.gdc.cancer.gov/files'
response = requests.get(files_endpt)
print(response.content)
```
