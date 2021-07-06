#create a rectangle latitude/longitude based mesh
#we will use some GIS way to define it
#longitude left and latitude bottom and nrow and ncolumn and resolution is used to define the rectangle
#because it is mesh, it represent the edge instead of center
#we will use gdal api for most operations
import os, sys
import numpy as np
from osgeo import ogr, osr, gdal, gdalconst
from shapely.geometry import Point, LineString, MultiLineString
from shapely.wkt import loads
from pystream.shared.hexagon import pyhexagon
from pystream.format.convert_coordinates_to_cell import convert_coordinates_to_cell

def create_hexagon_mesh(dX_left, dY_bot, dResolution_meter, ncolumn, nrow, sFilename_output, sFilename_shapefile):

    
    if os.path.exists(sFilename_output): 
        #delete it if it exists
        os.remove(sFilename_output)

    #pDriver = ogr.GetDriverByName('Esri Shapefile')
    pDriver = ogr.GetDriverByName('GeoJSON')
    #geojson
    pDataset = pDriver.CreateDataSource(sFilename_output)
    
    pDriver_shapefile = ogr.GetDriverByName('ESRI Shapefile')
    pDataset_shapefile = pDriver_shapefile.Open(sFilename_shapefile, 0)
    pLayer_shapefile = pDataset_shapefile.GetLayer(0)
    pSrs = pLayer_shapefile.GetSpatialRef()

    pLayer = pDataset.CreateLayer('cell', pSrs, ogr.wkbPolygon)
    # Add one attribute
    pLayer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger64)) #long type for high resolution
    
    pLayerDefn = pLayer.GetLayerDefn()
    pFeature = ogr.Feature(pLayerDefn)

    

    xleft = dX_left
    ybottom = dY_bot

    dArea = np.power(dResolution_meter,2.0)
    #hexagon edge
    dLength_edge = np.sqrt(  2.0 * dArea / (3.0* np.sqrt(3.0))  )
    dX_shift = 0.5 * dLength_edge
    dY_shift = 0.5 * dLength_edge * np.sqrt(3.0)
    dX_spacing = dLength_edge * 1.5
    dY_spacing = dLength_edge * np.sqrt(3.0)

    lID =0 

    #geojson
    aHexagon=list()
    #.........
    #(x2,y2)-----(x3,y3)
    #   |           |
    #(x1,y1)-----(x4,y4)
    #...............
    for column in range(0, ncolumn):
        for row in range(0, nrow):
            if column % 2 == 0 :
            #define a polygon here
                x1 = xleft + (column * dX_spacing)
                y1 = ybottom + (row * dY_spacing)
            else:
                x1 = xleft + (column * dX_spacing) #- dX_shift
                y1 = ybottom + (row * dY_spacing) - dY_shift


            x2 = x1 - dX_shift
            y2 = y1 + dY_shift

            x3 = x1 
            y3 = y1 + dY_shift * 2.0

            x4 = x1 + dLength_edge
            y4 = y1 + dY_shift * 2.0

            x5 = x4 + dX_shift
            y5 = y1 + dY_shift

            x6 = x1 + dLength_edge
            y6 = y1         
           

            ring = ogr.Geometry(ogr.wkbLinearRing)
            ring.AddPoint(x1, y1)
            ring.AddPoint(x2, y2)
            ring.AddPoint(x3, y3)
            ring.AddPoint(x4, y4)
            ring.AddPoint(x5, y5)
            ring.AddPoint(x6, y6)
            ring.AddPoint(x1, y1)
            pPolygon = ogr.Geometry(ogr.wkbPolygon)
            pPolygon.AddGeometry(ring)

            pFeature.SetGeometry(pPolygon)
            pFeature.SetField("id", lID)
            pLayer.CreateFeature(pFeature)

            lID = lID + 1

            dummy = loads( ring.ExportToWkt() )
            aCoords = dummy.exterior.coords
            dummy1= np.array(aCoords)
            pHexagon = convert_coordinates_to_cell(1, dummy1)
            aHexagon.append(pHexagon)

            pass
        
    pDataset = pLayer = pFeature  = None      



    return aHexagon



