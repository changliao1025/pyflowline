import os
from pathlib import Path
import json
import numpy as np
from osgeo import ogr, osr, gdal, gdalconst
import matplotlib.pyplot as plt
from matplotlib import cm
import cartopy.crs as ccrs
from shapely.wkt import loads
import matplotlib.path as mpath
import matplotlib.patches as mpatches
desired_proj = ccrs.Orthographic(central_longitude=-75, central_latitude=42, globe=None)
desired_proj = ccrs.PlateCarree()

def pyflowline_plot_flowline(oBasin_in, sVariable_in = None):
    sWorkspace_output_basin = oBasin_in.sWorkspace_output_basin 
    if sVariable_in is not None:
        if sVariable_in == 'flowline_filter_json':
            sFilename_json = oBasin_in.sFilename_flowline_filter_json
        else:
            if sVariable_in == 'flowline_simplified':
                sFilename_out = oBasin_in.sFilename_flowline_segment_index_before_intersect
                sFilename_json = os.path.join(sWorkspace_output_basin, sFilename_out)
            else:
                sFilename_out = oBasin_in.sFilename_flowline_final
                sFilename_json = os.path.join(sWorkspace_output_basin, sFilename_out)
            pass
    else:
        #default 
        sFilename_json = oBasin_in.sFilename_flowline_filter_json
    #convert existing flowline into the wgs83 system
    
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
            #patch = mpatches.PathPatch(path,  \
            #    lw=1, transform=ccrs.PlateCarree(), alpha=0.8, edgecolor= colours[lID])
            #ax.add_patch(patch)
            # plot control points and connecting lines
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