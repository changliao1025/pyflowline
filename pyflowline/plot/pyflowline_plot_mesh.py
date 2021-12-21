import sys, os, stat
import numpy as np
from pathlib import Path
import json

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
from matplotlib.collections import PatchCollection

import cartopy.crs as ccrs

desired_proj = ccrs.Orthographic(central_longitude=-75, central_latitude=42, globe=None)

#desired_proj = ccrs.PlateCarree()

def pyflowline_plot_mesh(oPyflowline_in):

    sWorkspace_output_case = oPyflowline_in.sWorkspace_output

    
    
    sFilename_json  =  oPyflowline_in.sFilename_mesh
    
    fig = plt.figure( dpi=300 )
    fig.set_figwidth( 12 )
    fig.set_figheight( 12 )
    ax = fig.add_axes([0.1, 0.15, 0.75, 0.6] , projection=desired_proj )
    
    ax.set_global()

    dLat_min = 90
    dLat_max = -90
    dLon_min = 180
    dLon_max = -180


    with open(sFilename_json) as json_file:
        data = json.load(json_file)        

        ncell = len(data)
        lID =0 
        for i in range(ncell):
            pcell = data[i]
            lCellID = int(pcell['lCellID'])
            lCellID_downslope = int(pcell['lCellID_downslope'])
            x_start=float(pcell['dLongitude_center_degree'])
            y_start=float(pcell['dLatitude_center_degree'])
            

            avertex = pcell['vVertex']
            nvertex = len(avertex)
            aLocation= np.full( (nvertex, 2), 0.0, dtype=float )
           
            
            for k in range(nvertex):
                aLocation[k,0] = avertex[k]['dLongitude_degree']
                aLocation[k,1] = avertex[k]['dLatitude_degree']

                if aLocation[k,0] > dLon_max:
                    dLon_max = aLocation[k,0]
                
                if aLocation[k,0] < dLon_min:
                    dLon_min = aLocation[k,0]
                
                if aLocation[k,1] > dLat_max:
                    dLat_max = aLocation[k,1]

                if aLocation[k,1] < dLat_min:
                    dLat_min = aLocation[k,1]

            
            polygon = mpatches.Polygon(aLocation, closed=True,  alpha=0.8, edgecolor=b,transform=ccrs.PlateCarree() )
            #aPatch.append(polygon)
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


