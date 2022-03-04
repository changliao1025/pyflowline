import os, sys

import numpy as np


from osgeo import ogr
from pyflowline.classes.tin import pytin

from pyflowline.formats.convert_coordinates import convert_pcs_coordinates_to_cell

def create_tin_mesh(dX_left_in, dY_bot_in, dResolution_meter_in, ncolumn_in, nrow_in,
sFilename_output_in, sFilename_spatial_reference_in):
     
    
    if os.path.exists(sFilename_output_in): 
        #delete it if it exists
        os.remove(sFilename_output_in)
    
    pDriver_shapefile = ogr.GetDriverByName('Esri Shapefile')

    
    pDataset = pDriver_shapefile.CreateDataSource(sFilename_output_in)
    
    pDataset_shapefile = pDriver_shapefile.Open(sFilename_spatial_reference_in, 0)
    pLayer_shapefile = pDataset_shapefile.GetLayer(0)
    pSrs = pLayer_shapefile.GetSpatialRef()

    #pSrs = osr.SpatialReference()  
    #pSrs.ImportFromEPSG(4326)    # WGS84 lat/lon

    pLayer = pDataset.CreateLayer('cell', pSrs, ogr.wkbPolygon)
    # Add one attribute
    pLayer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger64)) #long type for high resolution
    
    pLayerDefn = pLayer.GetLayerDefn()
    pFeature = ogr.Feature(pLayerDefn)

    

    xleft = dX_left_in
    ybottom = dY_bot_in

    dArea = np.power(dResolution_meter_in,2.0)
    #tin edge
    dLength_edge = np.sqrt(  4.0 * dArea /  np.sqrt(3.0) )  
    dX_shift = 0.5 * dLength_edge
    dY_shift = 0.5 * dLength_edge * np.sqrt(3.0) 
    dX_spacing = dX_shift * 2
    dY_spacing = dY_shift

    lID =0 

    #geojson
    aTin=list()
    #.........
    #(x2,y2)-----(x3,y3)
    #   |           |
    #(x1,y1)-----(x4,y4)
    #...............
    for column in range(0, ncolumn_in):
  
        for row in range(0, nrow_in):
            

            if column % 2 == 0 :
                if row % 2 == 0:
                    #define a polygon here
                    x1 = xleft + (column * dX_shift)
                    y1 = ybottom + (row * dY_spacing)
                    x2 = x1 + dX_spacing
                    y2 = y1 
                    x3 = x1 + dX_shift
                    y3 = y1 + dY_spacing
                else:
                    x1 = xleft + (column * dX_shift) 
                    y1 = ybottom + (row +1)* dY_spacing 
                    x2 = x1 + dX_shift
                    y2 = y1 - dY_shift    
                    x3 = x1 + dX_spacing
                    y3 = y1  
                    
            else:
                if row % 2 == 0:
                    x1 = xleft + column *  dX_shift
                    y1 = ybottom + (row + 1)* dY_spacing
                    x2 = x1 + dX_shift
                    y2 = y1 - dY_shift    
                    x3 = x1 + dX_spacing
                    y3 = y1   
                else:
                    x1 = xleft + column *  dX_shift
                    y1 = ybottom + (row )* dY_spacing
                    x2 = x1 + dX_spacing
                    y2 = y1   
                    x3 = x1 + dX_shift
                    y3 = y1 + dY_spacing
                         

                   
           
            aCoords = np.full((4,2), -9999.0, dtype=float)

            ring = ogr.Geometry(ogr.wkbLinearRing)
            ring.AddPoint(x1, y1)
            ring.AddPoint(x2, y2)
            ring.AddPoint(x3, y3)
            
            ring.AddPoint(x1, y1)
            pPolygon = ogr.Geometry(ogr.wkbPolygon)
            pPolygon.AddGeometry(ring)

            pFeature.SetGeometry(pPolygon)
            pFeature.SetField("id", lID)
            pLayer.CreateFeature(pFeature)

            lID = lID + 1

            #dummy = loads( ring.ExportToWkt() )
            #aCoords = dummy.exterior.coords
            aCoords[0,0] = x1
            aCoords[0,1] = y1
            aCoords[1,0] = x2
            aCoords[1,1] = y2
            aCoords[2,0] = x3
            aCoords[2,1] = y3
            
            aCoords[3,0] = x1
            aCoords[3,1] = y1
            
            dummy1= np.array(aCoords)
            pHexagon = convert_pcs_coordinates_to_cell(1, dummy1)
            aTin.append(pHexagon)

            pass
        
    pDataset = pLayer = pFeature  = None      


    return aTin
