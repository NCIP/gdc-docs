import requests
import json
import re

# This script will not work until $TOKEN_FILE_PATH  is replaced with an actual path.

with open("$TOKEN_FILE_PATH","r") as token:
    token_string = str(token.read().strip())

headers = {
           'X-Auth-Token': token_string
          }

data_endpt = 'https://api.gdc.cancer.gov/data/'
data_uuid = 'a1c1b23b-cc41-4e85-b1b7-62a42873c5af'
headers = {
           'X-Auth-Token': token_string
          }
response = requests.get(data_endpt + data_uuid, headers=headers)

# The file name can be found in the header within the Content-Disposition key.
response_head_cd = response.headers["Content-Disposition"]

file_name = re.findall("filename=(.+)", response_head_cd)[0]

with open(file_name, "wb") as output_file:
    output_file.write(response.content)