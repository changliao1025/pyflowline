#create a rectangle latitude/longitude based mesh
#we will use some GIS way to define it
#longitude left and latitude bottom and nrow and ncolumn and resolution is used to define the rectangle
#because it is mesh, it represent the edge instead of center
#we will use gdal api for most operations
import os, sys
import numpy as np
from osgeo import ogr, osr
from pyflowline.classes.latlon import pylatlon
from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_cell

def create_latlon_mesh(dLongitude_left_in, dLatitude_bot_in, dResolution_degree_in, ncolumn_in, nrow_in, sFilename_output_in):   
    if os.path.exists(sFilename_output_in): 
        os.remove(sFilename_output_in)

    #pDriver_shapefile = ogr.GetDriverByName('Esri Shapefile')
    pDriver_geojson = ogr.GetDriverByName('GeoJSON')    
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
    xleft = dLongitude_left_in
    xspacing= dResolution_degree_in
    ybottom = dLatitude_bot_in
    yspacing = dResolution_degree_in
    lCellID = 1
    #.........
    #(x2,y2)-----(x3,y3)
    #   |           |
    #(x1,y1)-----(x4,y4)
    #...............
    aLatlon = list()
    for iColumn in range(1, ncolumn_in+1):
        for iRow in range(1, nrow_in+1):
            #define a polygon here
            x1 = xleft + ((iColumn-1) * xspacing)
            y1 = ybottom + ((iRow-1) * yspacing)

            x2 = xleft + ((iColumn-1) * xspacing)
            y2 = ybottom + ((iRow ) * yspacing)

            x3 = xleft + ((iColumn ) * xspacing)
            y3 = ybottom + ((iRow ) * yspacing)

            x4 = xleft + ((iColumn ) * xspacing)
            y4 = ybottom + ((iRow-1) * yspacing)           

            ring = ogr.Geometry(ogr.wkbLinearRing)
            ring.AddPoint(x1, y1)
            ring.AddPoint(x2, y2)
            ring.AddPoint(x3, y3)
            ring.AddPoint(x4, y4)
            ring.AddPoint(x1, y1)
            pPolygon = ogr.Geometry(ogr.wkbPolygon)
            pPolygon.AddGeometry(ring)

            dLongitude_center = (x1 + x2 + x3 + x4)/4.0
            dLatitude_center = (y1 + y2 + y3 + y4)/4.0
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
            pLatlon = convert_gcs_coordinates_to_cell(3, dLongitude_center, dLatitude_center, dummy1)
            pLatlon.lCellID = lCellID
            dArea = pLatlon.calculate_cell_area()
            pLatlon.calculate_edge_length()
            pLatlon.dLongitude_center_degree = dLongitude_center
            pLatlon.dLatitude_center_degree = dLatitude_center
            pFeature.SetGeometry(pPolygon)
            pFeature.SetField("id", lCellID)
            pFeature.SetField("lon", dLongitude_center )
            pFeature.SetField("lat", dLatitude_center )
            pFeature.SetField("area", dArea )
            pLayer.CreateFeature(pFeature)

            lCellID_center = lCellID

            aNeighbor = list()
            if iRow > 1:#under
                lCellID0 = lCellID_center - 1
                aNeighbor.append(lCellID0)
                if iColumn > 1:
                    lCellID2 = lCellID0 - nrow_in
                    aNeighbor.append(lCellID2)

            if iColumn> 1:#left
                lCellID1 = nrow_in * (iColumn-2) + iRow 
                aNeighbor.append(lCellID1)  
                if iRow < nrow_in:
                    lCellID4 = lCellID1 + 1
                    aNeighbor.append(lCellID4)      
                    
            if iRow < nrow_in:#top
                lCellID3 = lCellID_center + 1
                aNeighbor.append(lCellID3)
                if iColumn < ncolumn_in:
                    lCellID6 = lCellID3 + nrow_in
                    aNeighbor.append(lCellID6) 
                    
            if iColumn  < ncolumn_in  : #right
                lCellID5 = nrow_in * iColumn + iRow 
                aNeighbor.append(lCellID5)
                if iRow > 1:
                    lCellID7 = lCellID5 -1
                    aNeighbor.append(lCellID7) 

            pLatlon.aNeighbor = aNeighbor
            pLatlon.nNeighbor = len(aNeighbor)
            pLatlon.aNeighbor_land= aNeighbor
            pLatlon.nNeighbor_land= pLatlon.nNeighbor            
            aLatlon.append(pLatlon)
            lCellID = lCellID + 1

            pass
        
    pDataset = pLayer = pFeature  = None      

    #calculate neighbor distance
    for pLatlon in aLatlon:
        aNeighbor = pLatlon.aNeighbor
        pLatlon.aNeighbor_distance=list()
        for lCellID1 in aNeighbor:
            for pLatlon1 in aLatlon:
                if pLatlon1.lCellID == lCellID1:
                    dDistance = pLatlon.pVertex_center.calculate_distance( pLatlon1.pVertex_center )
                    pLatlon.aNeighbor_distance.append(dDistance)
                    break

    return aLatlon





