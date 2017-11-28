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

print(response.content.decode("utf-8"))
