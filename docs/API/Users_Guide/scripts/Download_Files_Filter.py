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
