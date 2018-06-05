import json
import sys
import re

obo_file = sys.argv[1]
json_file = sys.argv[2]
txt = open(obo_file).read()
out = {}
for chunk in txt.split('[Term]'):
    id = re.search('id: (\w*:\w*)', chunk)
    name = re.search('name: (.*)', chunk)
    if id and name:
        out[id.group(1)] = name.group(1)
json.dump(out, open(json_file, 'w'))
