import requests
import json

token_file = <TOKEN_FILE_PATH>

file_id = "11443f3c-9b8b-4e47-b5b7-529468fec098"

data_endpt = 'https://api.gdc.cancer.gov/slicing/view/%s' % file_id

with open(token_file,'r') as token:
    token_string = str(token.read().strip())

params = { "gencode": ["BRCA1","BRCA2"] }

response = requests.post(data_endpt, data = json.dumps(params), headers = {"Content-Type": "application/json", "X-Auth-Token":token_string })

file_name = "brca_slices.bam"

with open(file_name,'wb') as output_file:
    output_file.write(response.content)
