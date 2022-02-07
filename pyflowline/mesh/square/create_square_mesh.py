#create a rectangle latitude/longitude based mesh
#we will use some GIS way to define it
#longitude left and latitude bottom and nrow and ncolumn and resolution is used to define the rectangle
#because it is mesh, it represent the edge instead of center
#we will use gdal api for most operations
import os, sys
from osgeo import ogr, osr
import numpy as np
from pyflowline.classes.square import pysquare
from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_cell


from pyflowline.algorithms.auxiliary.gdal_functions import  reproject_coordinates_batch

def create_square_mesh(dX_left_in, dY_bot_in, dResolution_meter_in, ncolumn_in, nrow_in, \
    sFilename_output_in, sFilename_spatial_reference_in):

   
    if os.path.exists(sFilename_output_in): 
        #delete it if it exists
        os.remove(sFilename_output_in)

    pDriver_shapefile = ogr.GetDriverByName('Esri Shapefile')
    pDriver_geojson = ogr.GetDriverByName('GeoJSON')

    pDataset_shapefile = pDriver_shapefile.Open(sFilename_spatial_reference_in, 0)
    pLayer_shapefile = pDataset_shapefile.GetLayer(0)
    pSpatial_reference = pLayer_shapefile.GetSpatialRef()   
        

    pDataset = pDriver_geojson.CreateDataSource(sFilename_output_in)
    
    pSpatial_reference_gcs = osr.SpatialReference()  
    pSpatial_reference_gcs.ImportFromEPSG(4326)    # WGS84 lat/lon     
    pSpatial_reference_gcs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

    pLayer = pDataset.CreateLayer('cell', pSpatial_reference_gcs, ogr.wkbPolygon)
    # Add one attribute
    pLayer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger64)) #long type for high resolution
    pLayer.CreateField(ogr.FieldDefn('lon', ogr.OFTReal)) #long type for high resolution
    pLayer.CreateField(ogr.FieldDefn('lat', ogr.OFTReal)) #long type for high resolution
    pArea_field = ogr.FieldDefn('area', ogr.OFTReal)
    pArea_field.SetWidth(20)
    pArea_field.SetPrecision(2)
    pLayer.CreateField(pArea_field)

    pLayerDefn = pLayer.GetLayerDefn()
    pFeature = ogr.Feature(pLayerDefn)    

    xleft = dX_left_in
    xspacing= dResolution_meter_in
    ybottom = dY_bot_in
    yspacing = dResolution_meter_in

    lID =0 
    #.........
    #(x2,y2)-----(x3,y3)
    #   |           |
    #(x1,y1)-----(x4,y4)
    #...............
    aSquare = list()
    for column in range(0, ncolumn_in):
        for row in range(0, nrow_in):
            #define a polygon here
            x1 = xleft + (column * xspacing)
            y1 = ybottom + (row * yspacing)

            x2 = xleft + (column * xspacing)
            y2 = ybottom + ((row + 1) * yspacing)

            x3 = xleft + ((column + 1) * xspacing)
            y3 = ybottom + ((row + 1) * yspacing)

            x4 = xleft + ((column + 1) * xspacing)
            y4 = ybottom + (row * yspacing)

           
            x = list()
            x.append(x1)
            x.append(x2)
            x.append(x3)
            x.append(x4)
          
            y = list()
            y.append(y1)
            y.append(y2)
            y.append(y3)
            y.append(y4)
           
            x_new , y_new = reproject_coordinates_batch(x, y, pSpatial_reference)
            x1=x_new[0]
            x2=x_new[1]
            x3=x_new[2]
            x4=x_new[3]
          
            y1=y_new[0]
            y2=y_new[1]
            y3=y_new[2]
            y4=y_new[3]        
           

            ring = ogr.Geometry(ogr.wkbLinearRing)
            ring.AddPoint(x1, y1)
            ring.AddPoint(x2, y2)
            ring.AddPoint(x3, y3)
            ring.AddPoint(x4, y4)
            ring.AddPoint(x1, y1)
            pPolygon = ogr.Geometry(ogr.wkbPolygon)
            pPolygon.AddGeometry(ring)

            dLon = (x1 + x2 + x3 + x4)/4.0
            dLat = (y1 + y2 + y3 + y4)/4.0
            aCoords = np.full((5,2), -9999.0, dtype=float)
            aCoords[0,0] = x1
            aCoords[0,1] = y1
            aCoords[1,0] = x2
            aCoords[1,1] = y2
            aCoords[2,0] = x3
            aCoords[2,1] = y3
            aCoords[3,0] = x4
            aCoords[3,1] = y4
            aCoords[4,0] = x1
            aCoords[4,1] = y1
            dummy1= np.array(aCoords)

            pSquare = convert_gcs_coordinates_to_cell(2, dLon, dLat, dummy1)
            dArea = pSquare.calculate_cell_area()

            pFeature.SetGeometry(pPolygon)
            pFeature.SetField("id", lID)
            pFeature.SetField("lon", dLon )
            pFeature.SetField("lat", dLat )
            pFeature.SetField("area", dArea )
            pLayer.CreateFeature(pFeature)

            lID = lID + 1
            aSquare.append(pSquare)


            pass

    pDataset = pLayer = pFeature  = None      



    return aSquare

