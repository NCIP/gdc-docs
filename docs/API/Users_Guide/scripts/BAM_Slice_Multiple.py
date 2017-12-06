import requests
import json

token_file = "$TOKEN_FILE_PATH"

file_ids = [
    "11443f3c-9b8b-4e47-b5b7-529468fec098",
    "1f103620-bb34-46f1-b565-94f0027e396d",
    "ca549554-a244-4209-9086-92add7bb7109"
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
