import numpy as np
import json
from json import JSONEncoder
class NumpyArrayEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return JSONEncoder.default(self, obj)

def export_mesh_to_json(aCell_in,sFilename_json_out):
    ncell=len(aCell_in)

    with open(sFilename_json_out, 'w', encoding='utf-8') as f:

        for i in range(1, ncell+1):
            sJson = aCell_in[i-1].export_to_json()
            #print(sJson)
            f.write(sJson)

        f.close()
    
    return
