#create a rectangle latitude/longitude based mesh
#we will use some GIS way to define it
#longitude left and latitude bottom and nrow and ncolumn and resolution is used to define the rectangle
#because it is mesh, it represent the edge instead of center
#we will use gdal api for most operations
import os, sys
from osgeo import ogr, osr, gdal, gdalconst



from pyflowline.algorithms.auxiliary.reproject_coordinates import reproject_coordinates, reproject_coordinates_batch

def create_square_mesh(dX_left_in, dY_bot_in, dResolution_meter_in, ncolumn_in, nrow_in, \
    sFilename_output_in, sFilename_spatial_reference_in):

   
    if os.path.exists(sFilename_output_in): 
        #delete it if it exists
        os.remove(sFilename_output_in)

    pDriver_shapefile = ogr.GetDriverByName('Esri Shapefile')
    #pDriver_geojson = ogr.GetDriverByName('GeoJSON')

    pDataset_shapefile = pDriver_shapefile.Open(sFilename_spatial_reference_in, 0)
    pLayer_shapefile = pDataset_shapefile.GetLayer(0)
    pSpatial_reference = pLayer_shapefile.GetSpatialRef()   
        

    pDataset = pDriver_shapefile.CreateDataSource(sFilename_output_in)
    
    pSpatial_reference_gcs = osr.SpatialReference()  
    pSpatial_reference_gcs.ImportFromEPSG(4326)    # WGS84 lat/lon     
    pSpatial_reference_gcs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

    pLayer = pDataset.CreateLayer('cell', pSpatial_reference_gcs, ogr.wkbPolygon)
    # Add one attribute
    pLayer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger64)) #long type for high resolution
    
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

            #x1,y1 = reproject_coordinates(x1, y1, pSpatial_reference)
            #x2,y2 = reproject_coordinates(x2, y2, pSpatial_reference)
            #x3,y3 = reproject_coordinates(x3, y3, pSpatial_reference)
            #x4,y4 = reproject_coordinates(x4, y4, pSpatial_reference)
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

            pFeature.SetGeometry(pPolygon)
            pFeature.SetField("id", lID)
            pLayer.CreateFeature(pFeature)

            lID = lID + 1


            pass
    pDataset = pLayer = pFeature  = None      



    return

