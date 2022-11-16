import json 
import requests
from dateutil import parser
from prose import Observatio

dates = []

for sha in versions:
    resp = requests.get(f"https://api.github.com/repos/lgrcia/prose/commits/{sha}")
    date = parser.parse(json.loads(resp.content)['commit']['author']['date'])

