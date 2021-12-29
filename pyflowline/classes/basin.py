import os
from abc import ABCMeta, abstractmethod
import json
from json import JSONEncoder
from pathlib import Path
import numpy as np
from osgeo import ogr, osr, gdal, gdalconst
from shapely.wkt import loads
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
from matplotlib import cm

import cartopy.crs as ccrs
desired_proj = ccrs.Orthographic(central_longitude=-75, central_latitude=42, globe=None)
desired_proj = ccrs.PlateCarree()

from pyflowline.formats.convert_shapefile_to_json import convert_shapefile_to_json

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
    sWorkspace_output_basin=''
    sFilename_flowline_raw=''    
    sFilename_flowline_filter=''
    sFilename_flowline_filter_json=''
    sFilename_dam=''
    sFilename_flowline_topo=''
    #before intersect
    sFilename_flowline_segment_order_before_intersect=''
    sFilename_flowline_segment_index_before_intersect=''
    sFilename_flowline_final=''
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
            self.sFilename_flowline_filter   = ''

        if 'sWorkspace_output_basin' in aParameter:
            self.sWorkspace_output_basin = aParameter['sWorkspace_output_basin']
        else:
            self.sWorkspace_output_basin   = ''
            print('The basin output path is not specified!')


        self.sFilename_flowline_filter_json = os.path.join(str(self.sWorkspace_output_basin ), "flowline_filter.json"  )

        if 'sFilename_dam' in aParameter:
            self.sFilename_dam = aParameter['sFilename_dam']
        else:
            self.sFilename_dam   = ''

        if 'sFilename_flowline_topo' in aParameter:
            self.sFilename_flowline_topo = aParameter['sFilename_flowline_topo']
        else:
            self.sFilename_flowline_topo   =''

        sBasinID = "{:03d}".format(self.lBasinID)

        self.sFilename_flowline_segment_index_before_intersect = 'flowline_segment_index_before_intersect_' + sBasinID + '.json'
        self.sFilename_flowline_segment_order_before_intersect = 'flowline_segment_order_before_intersect_' + sBasinID + '.json'
        self.sFilename_flowline_intersect  = 'flowline_intersect_' + sBasinID + '.json'
        self.sFilename_flowline_final = 'flowline_final_' + sBasinID + '.json'
        
        return
    
    def tojson(self):
        sJson = json.dumps(self.__dict__, \
            sort_keys=True, \
                indent = 4, \
                    ensure_ascii=True, \
                        cls=BasinClassEncoder)
        return sJson
    
    def convert_flowline_to_json(self):
        sFilename_raw = self.sFilename_flowline_filter            
        sFilename_out = self.sFilename_flowline_filter_json
        convert_shapefile_to_json(1, sFilename_raw, sFilename_out)
    
    def plot(self, sVariable_in=None):

        if sVariable_in is not None:
            if sVariable_in == 'flowline_filter_json':
                sFilename_json = self.sFilename_flowline_filter_json
            else:
                if sVariable_in == 'flowline_simplified':
                    sFilename_out = self.sFilename_flowline_segment_index_before_intersect
                    sFilename_json = os.path.join(sWorkspace_output_basin, sFilename_out)
                else:
                    sFilename_out = self.sFilename_flowline_final
                    sFilename_json = os.path.join(sWorkspace_output_basin, sFilename_out)
                pass
        else:
            #default 
            sFilename_json = self.sFilename_flowline_filter_json
        
        fig = plt.figure( dpi=150 )
        fig.set_figwidth( 4 )
        fig.set_figheight( 4 )
        ax = fig.add_axes([0.1, 0.15, 0.75, 0.8] , projection=desired_proj  )
        pDriver = ogr.GetDriverByName('GeoJSON')
        pDataset = pDriver.Open(sFilename_json, gdal.GA_ReadOnly)
        pLayer = pDataset.GetLayer(0)
    
        pSrs = osr.SpatialReference()  
        pSrs.ImportFromEPSG(4326)    # WGS84 lat/lon
    
        lID = 0
        dLat_min = 90
        dLat_max = -90
        dLon_min = 180
        dLon_max = -180
        n_colors = pLayer.GetFeatureCount()
        
        colours = cm.rainbow(np.linspace(0, 1, n_colors))
        for pFeature_shapefile in pLayer:
            pGeometry_in = pFeature_shapefile.GetGeometryRef()
            sGeometry_type = pGeometry_in.GetGeometryName()
            if sGeometry_type =='LINESTRING':
                dummy0 = loads( pGeometry_in.ExportToWkt() )
                aCoords_gcs = dummy0.coords
                aCoords_gcs= np.array(aCoords_gcs)
                nvertex = len(aCoords_gcs)
                for i in range(nvertex):
                    dLon = aCoords_gcs[i][0]
                    dLat = aCoords_gcs[i][1]
                    if dLon > dLon_max:
                        dLon_max = dLon
                    
                    if dLon < dLon_min:
                        dLon_min = dLon
                    
                    if dLat > dLat_max:
                        dLat_max = dLat
    
                    if dLat < dLat_min:
                        dLat_min = dLat
    
                codes = np.full(nvertex, mpath.Path.LINETO, dtype=int )
                codes[0] = mpath.Path.MOVETO
                path = mpath.Path(aCoords_gcs, codes)            
                x, y = zip(*path.vertices)
                line, = ax.plot(x, y, color= colours[lID])
                lID = lID + 1
                
    
        pDataset = pLayer = pFeature  = None      
    
        ax.set_extent([dLon_min  , dLon_max , dLat_min , dLat_max ])
        
        sDirname = os.path.dirname(sFilename_json)
        sFilename  = Path(sFilename_json).stem + '.png'
        sFilename_out = os.path.join(sDirname, sFilename)
        plt.savefig(sFilename_out, bbox_inches='tight')
        plt.show()
    
        return

    def preprocess_flowline(self):
        
        return