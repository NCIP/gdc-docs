import requests
import json

file_endpt = 'https://api.gdc.cancer.gov/files/'
file_uuid = 'd853e541-f16a-4345-9f00-88e03c2dc0bc'
response = requests.get(file_endpt + file_uuid)

# OUTPUT METHOD 1: Write to a file.
file = open("sample_request.json", "w")
file.write(response.text)
file.close()

# OUTPUT METHOD 2: View on screen.
print(json.dumps(response.json(), indent=2))