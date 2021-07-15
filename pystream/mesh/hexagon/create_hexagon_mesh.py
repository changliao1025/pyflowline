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

def checkIfDuplicates(listOfElems):
    ''' Check if given list contains any duplicates '''    
    iFlag_unique = 1
    for elem in listOfElems:
        if listOfElems.count(elem) > 1:
            iFlag_unique = 0
            break
        else:
            pass
    
    return iFlag_unique

def create_hexagon_mesh(iFlag_rotation, dX_left, dY_bot, dResolution_meter, ncolumn, nrow, sFilename_mesh_out, sFilename_spatial_reference_in):

    
    if os.path.exists(sFilename_mesh_out): 
        #delete it if it exists
        os.remove(sFilename_mesh_out)

    pDriver_shapefile = ogr.GetDriverByName('Esri Shapefile')
    #pDriver_geojson = ogr.GetDriverByName('GeoJSON')
    
    pDataset = pDriver_shapefile.CreateDataSource(sFilename_mesh_out)
    
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

    

    xleft = dX_left
    ybottom = dY_bot

    dArea = np.power(dResolution_meter,2.0)
    #hexagon edge
    dLength_edge = np.sqrt(  2.0 * dArea / (3.0* np.sqrt(3.0))  )
    

    

    #geojson
    aHexagon=list()
    #.........
    #(x2,y2)-----(x3,y3)
    #   |           |
    #(x1,y1)-----(x4,y4)
    #...............
  
    lCellID = 1
    if iFlag_rotation ==0:
        dX_shift = 0.5 * dLength_edge * np.sqrt(3.0)
        dY_shift = 0.5 * dLength_edge
        dX_spacing = dLength_edge * np.sqrt(3.0)
        dY_spacing = dLength_edge * 1.5
        for iRow in range(1, nrow+1):
            for iColumn in range(1, ncolumn+1):
                if iRow % 2 == 1 : #odd
                #define a polygon here
                    x1 = xleft + (iColumn-1) * dX_spacing
                    y1 = ybottom + (iRow-1) * dY_spacing
                else:
                    x1 = xleft + (iColumn-1) * dX_spacing + dX_shift
                    y1 = ybottom + (iRow-1) * dY_spacing
    
    
                x2 = x1 
                y2 = y1 + dLength_edge
    
                x3 = x1 + dX_shift
                y3 = y2 + dY_shift
    
                x4 = x1 + dX_spacing
                y4 = y2
    
                x5 = x4 
                y5 = y1 
    
                x6 = x3
                y6 = y1  -   dY_shift     
               
                aCoords = np.full((7,2), -9999.0, dtype=float)
    
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
                pFeature.SetField("id", lCellID)
                pLayer.CreateFeature(pFeature)
                
    
    
                #dummy = loads( ring.ExportToWkt() )
                #aCoords = dummy.exterior.coords
                aCoords[0,0] = x1
                aCoords[0,1] = y1
                aCoords[1,0] = x2
                aCoords[1,1] = y2
                aCoords[2,0] = x3
                aCoords[2,1] = y3
                aCoords[3,0] = x4
                aCoords[3,1] = y4
                aCoords[4,0] = x5
                aCoords[4,1] = y5
                aCoords[5,0] = x6
                aCoords[5,1] = y6
                aCoords[6,0] = x1
                aCoords[6,1] = y1
           
                dummy1= np.array(aCoords)
                pHexagon = convert_coordinates_to_cell(1, dummy1)
                pHexagon.lCellID = lCellID
                lCellID_center = lCellID
                #build topoloy
                aNeighbor=list()
                if iColumn > 1:#0
                    lCellID0 = lCellID_center - 1
                    aNeighbor.append(lCellID0)

                if iRow < nrow :#1 and 2
                    if iRow %2 ==0:
                        lCellID1 = ncolumn * iRow + iColumn 
                        aNeighbor.append(lCellID1)
                        if iColumn!=ncolumn:
                            lCellID2 = ncolumn * iRow + iColumn + 1
                            aNeighbor.append(lCellID2)
                    else:
                        lCellID2 = ncolumn * iRow + iColumn
                        aNeighbor.append(lCellID2)
                        if iColumn != 1:
                            lCellID1 = ncolumn * iRow + iColumn - 1 
                            aNeighbor.append(lCellID1)

                
                        
                        
                if iColumn < ncolumn:#3
                    lCellID3 = lCellID_center + 1
                    aNeighbor.append(lCellID3)

                if iRow > 1 : #4 and 5
                    if iRow %2 ==1:
                        lCellID4 = ncolumn * (iRow-2) + iColumn 
                        aNeighbor.append(lCellID4)

                        if iColumn !=1:
                            lCellID5 = ncolumn * (iRow-2) + iColumn -1
                            aNeighbor.append(lCellID5)
                    else:
                        lCellID5 = ncolumn * (iRow-2) + iColumn 
                        aNeighbor.append(lCellID5)
                        if iColumn!=ncolumn:
                            lCellID4 = ncolumn * (iRow-2) + iColumn + 1
                            aNeighbor.append(lCellID4)


                if checkIfDuplicates(aNeighbor) == 0:
                    print('error')        


                pHexagon.aNeighbor = aNeighbor
                print(lCellID_center,aNeighbor)
                pHexagon.nNeighbor = len(aNeighbor)
                aHexagon.append(pHexagon)
                lCellID= lCellID + 1    
                pass
              
        pass
    else:
        dX_shift = 0.5 * dLength_edge
        dY_shift = 0.5 * dLength_edge * np.sqrt(3.0)
        dX_spacing = dLength_edge * 1.5
        dY_spacing = dLength_edge * np.sqrt(3.0)
        for iColumn in range(1, ncolumn+1):
            for iRow in range(1, nrow+1):
                if iColumn % 2 == 0 :
                #define a polygon here
                    x1 = xleft + (iColumn-1) * dX_spacing
                    y1 = ybottom + (iRow-1) * dY_spacing
                else:
                    x1 = xleft + (iColumn-1) * dX_spacing #- dX_shift
                    y1 = ybottom + (iRow-1) * dY_spacing + dY_shift
    
    
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
               
                aCoords = np.full((7,2), -9999.0, dtype=float)
    
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
                pFeature.SetField("id", lCellID)
                pLayer.CreateFeature(pFeature)
                
    
    
                #dummy = loads( ring.ExportToWkt() )
                #aCoords = dummy.exterior.coords
                aCoords[0,0] = x1
                aCoords[0,1] = y1
                aCoords[1,0] = x2
                aCoords[1,1] = y2
                aCoords[2,0] = x3
                aCoords[2,1] = y3
                aCoords[3,0] = x4
                aCoords[3,1] = y4
                aCoords[4,0] = x5
                aCoords[4,1] = y5
                aCoords[5,0] = x6
                aCoords[5,1] = y6
                aCoords[6,0] = x1
                aCoords[6,1] = y1
                    
                dummy1= np.array(aCoords)
                pHexagon = convert_coordinates_to_cell(1, dummy1)
                pHexagon.lCellID = lCellID
                #build topoloy
                lCellID_center = lCellID
                
                aNeighbor=list()
                if iRow > 1:#0
                    lCellID0 = lCellID_center - 1
                    aNeighbor.append(lCellID0)

                if iColumn> 1:#1 ans 2
                    if iColumn %2 ==0:
                        lCellID1 = nrow * (iColumn-2) + iRow 
                        aNeighbor.append(lCellID1)
                        if iRow!=nrow:
                            lCellID2 = nrow * (iColumn-2) + iRow +1
                            aNeighbor.append(lCellID2)
                    else:
                        lCellID2 = nrow * (iColumn-2) + iRow
                        aNeighbor.append(lCellID2)
                        if iRow != 1:
                            lCellID1 = nrow * (iColumn-2) + iRow -1
                            aNeighbor.append(lCellID1)

                                        
                        
                if iRow < nrow:#3
                    lCellID3 = lCellID_center + 1
                    aNeighbor.append(lCellID3)

                if iColumn  < ncolumn  : #4 and 5
                    if iColumn %2 ==1:
                        lCellID4 = nrow * iColumn + iRow 
                        aNeighbor.append(lCellID4)
                        if iRow !=1:
                            lCellID5 = nrow * iColumn + iRow -1
                            aNeighbor.append(lCellID5)
                    else:
                        lCellID5 = nrow * iColumn + iRow 
                        aNeighbor.append(lCellID5)
                        if iRow!=nrow:
                            lCellID4 = nrow * iColumn + iRow  +1
                            aNeighbor.append(lCellID4)
                                                   
                if checkIfDuplicates(aNeighbor) == 0:
                    print('error')   
                pHexagon.aNeighbor = aNeighbor
                print(lCellID_center,aNeighbor)
                pHexagon.nNeighbor = len(aNeighbor)
                aHexagon.append(pHexagon)
                lCellID= lCellID +1
    
                pass
       

        
        
    pDataset = pLayer = pFeature  = None      
    

    



    return aHexagon



