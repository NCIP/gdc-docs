import requests
import json
import re

data_endpt = "https://api.gdc.cancer.gov/data"

ids = [
    "b658d635-258a-4f6f-8377-767a43771fe4",
    "3968213d-b293-4b3d-8033-5b5a0ca07b6c"
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
