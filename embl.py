import requests, sys

requestURL = "https://www.ebi.ac.uk/proteins/api/variation?offset=0&size=100&disease=alzheimer"

r = requests.get(requestURL, headers={ "Accept" : "text/x-gff"})

if not r.ok:
  r.raise_for_status()
  sys.exit()

responseBody = r.text
print(responseBody)
