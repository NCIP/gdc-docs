import requests
import json
import re

data_endpt = "https://api.gdc.cancer.gov/data"

ids = [
    "a68e75c8-0e01-446b-86b8-09b6903fa409",
    "649e1e0e-b4ce-491c-94c3-64273680159b"
    ]

params = {"ids": ids}

response = requests.post(data_endpt,
                        data = json.dumps(params),
                        headers={
                            "Content-Type": "application/json"
                            })

response_head_cd = response.headers["Content-Disposition"]

file_name = re.findall("filename=(.+)", response_head_cd)[0]

with open(file_name, "wb") as output_file:
    output_file.write(response.content)
