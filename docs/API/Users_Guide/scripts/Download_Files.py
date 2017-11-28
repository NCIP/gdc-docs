import requests
import json

cases_endpt = 'https://api.gdc.cancer.gov/data'

fields = ["submitter_id","case_id","primary_site","disease_type","diagnoses.vital_status"]
fields = ','.join(fields)

params = {
    "fields":fields,
    "format":"TSV",
    "size":"100",
    }

response = requests.get(cases_endpt, params = params)

print(response.content)
