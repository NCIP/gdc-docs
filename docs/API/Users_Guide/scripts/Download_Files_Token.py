import requests
import json
import re

token_file = <TOKEN_FILE_PATH>

file_id = "2f97081c-7e84-4a93-91a8-fee860769f8e"

data_endpt = 'https://api.gdc.cancer.gov/data/%s' % file_id

with open(token_file,'r') as token:
    token_string = str(token.read().strip())

response = requests.get(data_endpt, headers = {"Content-Type": "application/json", "X-Auth-Token":token_string })

response_head = response.headers["Content-Disposition"]

file_name = re.findall("filename=(.+)", response_head)[0]

with open(file_name,'wb') as output_file:
    output_file.write(response.content)
