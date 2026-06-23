import requests
status_endpt = "https://api.gdc.cancer.gov/status"
response = requests.get(status_endpt)

# OUTPUT METHOD 1: Write to a file.
file = open("api_status.json", "w")
file.write(response.text)
file.close()

# OUTPUT METHOD 2: View on screen.
print(response.content)