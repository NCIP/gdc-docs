import requests
import json
import re

'''
 This script will not work until $TOKEN_FILE_PATH
 is replaced with an actual path.
'''
token_file = "$TOKEN_FILE_PATH"

file_id = "2f97081c-7e84-4a93-91a8-fee860769f8e"

data_endpt = "https://api.gdc.cancer.gov/data/{}".format(file_id)

with open(token_file, "r") as token:
    token_string = str(token.read().strip())

response = requests.get(data_endpt,
                        headers = {
                            "Content-Type": "application/json",
                            "X-Auth-Token": token_string
                            })

response_head_cd = response.headers["Content-Disposition"]

file_name = re.findall("filename=(.+)", response_head_cd)[0]

with open(file_name, "wb") as output_file:
    output_file.write(response.content)
