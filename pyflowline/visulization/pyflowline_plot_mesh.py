import sys, os, stat
import numpy as np
from pathlib import Path
from osgeo import ogr, osr, gdal, gdalconst
import json
from shapely.wkt import loads
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
from matplotlib.collections import PatchCollection

import cartopy.crs as ccrs

desired_proj = ccrs.Orthographic(central_longitude=-75, central_latitude=42, globe=None)

desired_proj = ccrs.PlateCarree()

def pyflowline_plot_mesh(oPyflowline_in):

    sWorkspace_output_case = oPyflowline_in.sWorkspace_output
    
    sFilename_json  =  oPyflowline_in.sFilename_mesh
    
    fig = plt.figure( dpi=150 )
    fig.set_figwidth( 4 )
    fig.set_figheight( 4 )
    ax = fig.add_axes([0.1, 0.15, 0.75, 0.7] , projection=desired_proj )
    
    ax.set_global()
    pDriver = ogr.GetDriverByName('GeoJSON')
    pDataset = pDriver.Open(sFilename_json, gdal.GA_ReadOnly)
    pLayer = pDataset.GetLayer(0)
   
    pSrs = osr.SpatialReference()  
    pSrs.ImportFromEPSG(4326)    # WGS84 lat/lon

    dLat_min = 90
    dLat_max = -90
    dLon_min = 180
    dLon_max = -180


    for pFeature_shapefile in pLayer:
        pGeometry_in = pFeature_shapefile.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()

        
        lID =0 
        if sGeometry_type =='POLYGON':
            dummy0 = loads( pGeometry_in.ExportToWkt() )
            aCoords_gcs = dummy0.exterior.coords
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

            
            polygon = mpatches.Polygon(aCoords_gcs[:,0:2], closed=True,  \
                alpha=0.8, edgecolor = 'black',facecolor='none', \
                    transform=ccrs.PlateCarree() )
            
            ax.add_patch(polygon)                   
                    
    
    dDiff_lon = dLon_max - dLon_min
    dDiff_lat = dLat_max - dLat_min
   
    ax.set_extent([dLon_min  , dLon_max , dLat_min , dLat_max ])

    ax.coastlines()#resolution='110m')
    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.3, linestyle='--')

    
    sDirname = os.path.dirname(sFilename_json)
    sFilename  = Path(sFilename_json).stem + '.png'
    sFilename_out = os.path.join(sDirname, sFilename)
    plt.savefig(sFilename_out, bbox_inches='tight')
    
    pDataset = pLayer = pFeature  = None   
    plt.show()   
    return


