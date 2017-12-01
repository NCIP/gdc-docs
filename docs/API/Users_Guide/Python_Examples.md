# Using Python to Query the GDC API

Python can be a versatile tool for retrieving information from the GDC API and performing downstream processing. This guide details some examples that demonstrate the basic methods for querying the API using Python. The examples in this guide will use the [requests](http://docs.python-requests.org/en/master/) Python library and should be compatible with Python3.

## Querying Metadata

### A Basic Query

This example passes some basic parameters (fields, format, size) to the `cases` endpoint and prints the results. Note that the `fields` parameter needs to be a string comprising comma-delimited field names.  

```TXT
Choose the Python tab to view script.
```
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
[Download Script](scripts/Basic_Query.py)

### A Filtered Query

In the next example, a `filters` parameter is added to the query. This parameter is passed as a Python dictionary object. The filter used in this example will only show cases that come from a patient with kidney disease (`primary_site: Kidney`).

```TXT
Choose the Python tab to view script.
```
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

# The filters parameter needs to be converted from a dictionary to JSON-formatted string

params = {
    "filters":json.dumps(filters),
    "fields":fields,
    "format":"TSV",
    "size":"100",
    }

response = requests.get(cases_endpt, params = params)

print(response.content)
```
[Download Script](scripts/Filter_Query.py)
### Complex Filters

The following example utilizes the `and` operator in the filter to returns information about RNA-Seq BAM files that originate from lung cancer patients. Note that these three filters are nested within a list in the highest level `content` key.  

```TXT
Choose the Python tab to view script.
```
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
[Download Script](scripts/Complex_Query.py)


## Downloading Files

GDC files can also be downloaded from the API using Python scripts.

### A Simple Download Request

Here an open-access file is downloaded from the GDC using the file UUID.  

```TXT
Choose the Python tab to view script.
```
```Python
import requests
import json
import re

file_id = "b658d635-258a-4f6f-8377-767a43771fe4"

data_endpt = 'https://api.gdc.cancer.gov/data/%s' % file_id

response = requests.get(data_endpt, headers={"Content-Type": "application/json"})

response_head = response.headers["Content-Disposition"]

file_name = re.findall("filename=(.+)", response_head)[0]

with open(file_name,'wb') as output_file:
    output_file.write(response.content)
```
[Download Script](scripts/Download_Files_Post.py)

### Passing a Token to Download a Controlled-Access File

Here a token is passed to the script by specifying a plain text file that contains only the GDC token.  A token can be downloaded by logging into the GDC Data Portal. See the [Data Security](../../Data/Data_Security.md) documentation for more details.  

```TXT
Choose the Python tab to view script.
```
```Python
import requests
import json
import re

# The GDC token is passed by reading in the text file containing the token.
token_file = <TOKEN_FILE_PATH>

file_id = "2f97081c-7e84-4a93-91a8-fee860769f8e"

data_endpt = 'https://api.gdc.cancer.gov/data/%s' % file_id

# The token is the read into the token_string variable
with open(token_file,'r') as token:
    token_string = str(token.read().strip())

# The location of the token needs to be specified in the headers dictionary
response = requests.get(data_endpt, headers={"Content-Type": "application/json", "X-Auth-Token":token_string })

# The file name is located in the headers of the response
response_head = response.headers["Content-Disposition"]
file_name = re.findall("filename=(.+)", response_head)[0]

# Note that 'wb' is required for the script to work with Python3
with open(file_name,'wb') as output_file:
    output_file.write(response.content)
```
[Download Script](scripts/Download_Files_Token.py)

### Post Request to Download Multiple Files  

This example uses a Python list to specify a set of file UUIDs.  Note that the list in the example was populated manually but could potentially be populated from an external list.  

```TXT
Choose the Python tab to view script.
```
```Python
import requests
import json
import re

data_endpt = 'https://api.gdc.cancer.gov/data'

ids = ["b658d635-258a-4f6f-8377-767a43771fe4","3968213d-b293-4b3d-8033-5b5a0ca07b6c"]

params = {
    "ids":ids
    }

response = requests.post(data_endpt, data = json.dumps(params), headers={"Content-Type": "application/json"})

response_head = response.headers["Content-Disposition"]

file_name = re.findall("filename=(.+)", response_head)[0]

with open(file_name,'wb') as output_file:
    output_file.write(response.content)
```
[Download Script](scripts/Download_Files_Post.py)

### Downloading a Set of Files Based on a Filter

Here files based on a set of filters are downloaded.  First file UUIDs are retrieved based on the filters.  These UUIDs are required to download the correct files.   

```TXT
Choose the Python tab to view script.
```
```Python
import requests
import json
import re

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
            "field":"cases.demographic.race",
            "value":["white"]
            }
        },
        {
        "op":"in",
        "content":{
            "field":"cases.demographic.gender",
            "value":["female"]
            }
        },
        {
        "op":"in",
        "content":{
            "field":"files.analysis.workflow_type",
            "value":["HTSeq - FPKM"]
            }
        }
    ]
}

params = {
    "filters":json.dumps(filters),
    "fields":'file_id',
    "format":"JSON",
    "size":"1000"
    }

response = requests.get(files_endpt, params = params)

file_uuid_list = []

for file_entry in json.loads(response.content.decode("utf-8"))["data"]["hits"]:
    file_uuid_list.append(file_entry["file_id"])

data_endpt = 'https://api.gdc.cancer.gov/data'

params = {
    "ids":file_uuid_list
    }

response = requests.post(data_endpt, data = json.dumps(params), headers={"Content-Type": "application/json"})

response_head = response.headers["Content-Disposition"]

file_name = re.findall("filename=(.+)", response_head)[0]

with open(file_name,'wb') as output_file:
    output_file.write(response.content)
```
[Download Script](scripts/Download_Files_Filter.py)

## Querying and Parsing Data

The gene, mutation, and survival data used to generate the visualization data on the GDC Data Portal can also be retrieved from the API.  Below are instructions for doing this using Python.  





## Basic Troubleshooting

The following script should produce an unformatted JSON string when run. Run this to verify that a valid connection is being made to the GDC API.  

```Python
import requests
files_endpt = 'https://api.gdc.cancer.gov/files'
response = requests.get(files_endpt)
print(response.content)
```
