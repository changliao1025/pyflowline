
import json

from pyflowline.shared.basin import pybasin

aBasin = list()

nbasin = 13
aParameter={}
for i in range(nbasin):
    aParameter['lBasinID'] = i + 1
    pBasin = pybasin(aParameter)
    aBasin.append(pBasin)

sFilename_json_out = '/qfs/people/liao313/workspace/python/pyflowline/pyflowline/config/pyflowline_icom_basin.json'
with open(sFilename_json_out, 'w', encoding='utf-8') as f:
    sJson = json.dumps([json.loads(ob.tojson()) for ob in aBasin],\
         sort_keys=True, \
         indent = 4)        
    f.write(sJson)    
    f.close()