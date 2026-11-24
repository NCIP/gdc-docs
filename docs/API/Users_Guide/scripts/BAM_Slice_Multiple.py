import requests
import json

'''
 This script will not work until $TOKEN_FILE_PATH
 is replaced with an actual path.
'''
token_file = "$TOKEN_FILE_PATH"

file_ids = [
    "a5b9acbd-adea-47d0-a3a3-3eb6ebfde56b",
    "8b06714b-a39c-4a24-897d-170af8837b57",
    "4400e648-8c19-4707-9665-dbeb907b166c"
    ]

for file_id in file_ids:

    data_endpt = "https://api.gdc.cancer.gov/slicing/view/{}".format(file_id)

    with open(token_file, "r") as token:
        token_string = str(token.read().strip())

    params = {
        "regions": ["chr1:1-20000", "chr10:129000-160000"]
        }

    response = requests.post(data_endpt,
                            data = json.dumps(params),
                            headers = {
                                "Content-Type": "application/json",
                                "X-Auth-Token": token_string
                                })

    file_name = "{}_region_slices.bam".format(file_id)

    with open(file_name, "wb") as output_file:
        output_file.write(response.content)
