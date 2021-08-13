import numpy as np
import json

def export_mesh_info_to_json(aCell_in,sFilename_json_out):
    ncell=len(aCell_in)

    with open(sFilename_json_out, 'w', encoding='utf-8') as f:
        sJson = json.dumps([json.loads(ob.tojson()) for ob in aCell_in], indent = 4)        
        f.write(sJson)    
        f.close()
    
    return
