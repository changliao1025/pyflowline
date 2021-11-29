from abc import ABCMeta, abstractmethod
import numpy as np
import json
from json import JSONEncoder

class BasinClassEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.float):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        
        return JSONEncoder.default(self, obj)


class pybasin(object):
    lBasinID =1 
    lCellID_outlet=-1
    iFlag_disconnected =0
    iFlag_dam=0
    dLongitude_outlet_degree = -9999.
    dLatitude_outlet_degree = -9999.

    dAccumulation_threshold= 100000.0

    dThreshold_small_river = 10000

    sFilename_flowline_raw=''
    sFilename_flowline_filter=''
    sFilename_dam=''
    sFilename_flowline_topo=''
    #before intersect
    sFilename_flowline_segment_order_before_intersect=''
    sFilename_flowline_segment_index_before_intersect=''
    def __init__(self, aParameter):

        if 'lBasinID' in aParameter:            
            self.lBasinID             = int(aParameter['lBasinID'])
        else:
            self.lBasinID   = 1
        
        if 'lCellID_outlet' in aParameter:            
            self.lCellID_outlet             = int(aParameter['lCellID_outlet'])
        else:
            self.lCellID_outlet   = -1

        if 'iFlag_disconnected' in aParameter:            
            self.iFlag_disconnected             = int(aParameter['iFlag_disconnected'])
        else:
            self.iFlag_disconnected   = 0
        
        if 'iFlag_dam' in aParameter:            
            self.iFlag_dam             = int(aParameter['iFlag_dam'])
        else:
            self.iFlag_dam   = 0
        
        if 'dLongitude_outlet_degree' in aParameter:            
            self.dLongitude_outlet_degree             = float(aParameter['dLongitude_outlet_degree'])
        else:
            self.dLongitude_outlet_degree   = -9999.
        
        if 'dLatitude_outlet_degree' in aParameter:            
            self.dLatitude_outlet_degree             = float(aParameter['dLatitude_outlet_degree'])
        else:
            self.dLatitude_outlet_degree   = -9999.
        
        if 'dThreshold_small_river' in aParameter:            
            self.dThreshold_small_river             = float(aParameter['dThreshold_small_river'])
        else:
            self.dThreshold_small_river   =10000.0

        if 'dAccumulation_threshold' in aParameter:            
            self.dAccumulation_threshold             = float(aParameter['dAccumulation_threshold'])
        else:
            self.dAccumulation_threshold = 100000.0

        if 'sFilename_flowline_raw' in aParameter:
            self.sFilename_flowline_raw = aParameter['sFilename_flowline_raw']
        else:
            self.sFilename_flowline_raw   =''
       
        if 'sFilename_flowline_filter' in aParameter:
            self.sFilename_flowline_filter = aParameter['sFilename_flowline_filter']
        else:
            self.sFilename_flowline_filter   =''

        if 'sFilename_dam' in aParameter:
            self.sFilename_dam = aParameter['sFilename_dam']
        else:
            self.sFilename_dam   = ''

        if 'sFilename_flowline_topo' in aParameter:
            self.sFilename_flowline_topo = aParameter['sFilename_flowline_topo']
        else:
            self.sFilename_flowline_topo   =''

        sBasinID = "{:03d}".format(self.lBasinID)

        self.sFilename_flowline_segment_index_before_intersect = 'flowline_segment_index_before_intersect_' + sBasinID + '.shp'
        self.sFilename_flowline_segment_order_before_intersect = 'flowline_segment_order_before_intersect_' + sBasinID + '.shp'
        self.sFilename_flowline_intersect  = 'flowline_intersect_' + sBasinID + '.shp'
        
        return
    
    def tojson(self):
        sJson = json.dumps(self.__dict__, \
            sort_keys=True, \
                indent = 4, \
                    ensure_ascii=True, \
                        cls=BasinClassEncoder)
        return sJson