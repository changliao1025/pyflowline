#create a rectangle latitude/longitude based mesh
#we will use some GIS way to define it
#longitude left and latitude bottom and nrow and ncolumn and resolution is used to define the rectangle
#because it is mesh, it represent the edge instead of center
#we will use gdal api for most operations
import os
from osgeo import ogr, osr
import numpy as np
from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_cell
from pyflowline.external.pyearth.gis.gdal.gdal_functions import  reproject_coordinates_batch

def index_to_row_col(index, num_columns):
    index -= 1  # Adjust for 1-based indexing
    row = index // num_columns + 1
    col = index % num_columns + 1
    return row, col

def create_square_mesh(dX_left_in, dY_bot_in,
                        dResolution_meter_in,
                        ncolumn_in, nrow_in,
    sFilename_output_in, 
    sFilename_spatial_reference_in, 
    pBoundary_in):   
    """
    _summary_

    Args:
        dX_left_in (_type_): _description_
        dY_bot_in (_type_): _description_
        dResolution_meter_in (_type_): _description_
        ncolumn_in (_type_): _description_
        nrow_in (_type_): _description_
        sFilename_output_in (_type_): _description_
        sFilename_spatial_reference_in (_type_): _description_

    Returns:
        _type_: _description_
    """
    #for the reason that a geometry object will be crash if the associated dataset is closed, we must pass wkt string
    #https://gdal.org/api/python_gotchas.html
    pBoundary = ogr.CreateGeometryFromWkt(pBoundary_in)

    if os.path.exists(sFilename_output_in): 
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
    pLayer.CreateField(ogr.FieldDefn('cellid', ogr.OFTInteger64)) #long type for high resolution
    pLayer.CreateField(ogr.FieldDefn('longitude', ogr.OFTReal)) #long type for high resolution
    pLayer.CreateField(ogr.FieldDefn('latitude', ogr.OFTReal)) #long type for high resolution
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

    def add_cell_into_list(aList, lCellID, iRow, iColumn, dLongitude_center, dLatitude_center, aCoords ):          
    
        pSquare = convert_gcs_coordinates_to_cell(2, dLongitude_center, dLatitude_center, aCoords)
        pSquare.lCellID = lCellID
        dArea = pSquare.calculate_cell_area()
        pSquare.calculate_edge_length()                
        #build topoloy
        aNeighbor=list()
        aNeighbor_distance=list()
        #lCellID_center = lCellID
        #counter-clock wise direction to add the neighbor
        if iRow > 1:#under
            iRow_dummy = iRow - 1
            if iColumn > 1:
                iColumn_dummy = iColumn - 1
                lCellID2 = (iRow_dummy-1) * ncolumn_in + iColumn_dummy #lCellID0 - nrow_in
                aNeighbor.append(lCellID2)
                
            lCellID0 =  (iRow_dummy-1) * ncolumn_in + iColumn
            aNeighbor.append(lCellID0)                    
        if iColumn  < ncolumn_in  : #right
            iColumn_dummy = iColumn + 1
            if iRow > 1:
                iRow_dummy = iRow - 1
                lCellID7 = (iRow_dummy-1) * ncolumn_in + iColumn_dummy# lCellID5 -1
                aNeighbor.append(lCellID7) 
            lCellID5 = (iRow-1) * ncolumn_in + iColumn_dummy #nrow_in * iColumn + iRow 
            aNeighbor.append(lCellID5)                       
        if iRow < nrow_in:#top
            iRow_dummy = iRow + 1
            if iColumn < ncolumn_in:
                iColumn_dummy = iColumn + 1
                lCellID6 = (iRow_dummy-1) * ncolumn_in + iColumn_dummy #lCellID3 + nrow_in
                aNeighbor.append(lCellID6) 
            lCellID3 = (iRow_dummy-1) * ncolumn_in + iColumn #lCellID_center + 1
            aNeighbor.append(lCellID3)         
        
        if iColumn> 1:#left
            iColumn_dummy = iColumn - 1
            if iRow < nrow_in:
                iRow_dummy = iRow + 1
                lCellID4 = (iRow_dummy-1) * ncolumn_in + iColumn_dummy #lCellID1 + 1
                aNeighbor.append(lCellID4)   
            lCellID1 = (iRow-1) * ncolumn_in + iColumn_dummy #nrow_in * (iColumn-2) + iRow 
            aNeighbor.append(lCellID1)  

        pSquare.aNeighbor = aNeighbor
        pSquare.nNeighbor = len(aNeighbor)
        pSquare.aNeighbor_land= aNeighbor
        pSquare.nNeighbor_land= pSquare.nNeighbor
        aList.append(pSquare)

        return aList, dArea

    #change the order because mpas uses counter-clock wise to store the vertices
    #we will also start from the lower-left corner, and then go to the right and then go up
    #so the final index will be like this
    #3 4
    #1 2
    #lCellID = 1
    #.........
    #(x4,y4)-----(x3,y3)
    #   |           |
    #(x1,y1)-----(x2,y2)
    #...............
    aSquare = list()
    for iRow in range(1, nrow_in+1):    
        for iColumn in range(1, ncolumn_in+1):        
            #global cell id for the mesh
            lCellID = (iRow-1) * ncolumn_in + iColumn

            #define a polygon here
            x1 = xleft + ((iColumn-1) * xspacing)
            y1 = ybottom + ((iRow-1) * yspacing)

            x2 = xleft + ((iColumn ) * xspacing)
            y2 = ybottom + ((iRow-1) * yspacing)     

            x3 = xleft + ((iColumn ) * xspacing)
            y3 = ybottom + ((iRow ) * yspacing)

            x4 = xleft + ((iColumn-1) * xspacing)
            y4 = ybottom + ((iRow ) * yspacing)  

            x = [x1, x2, x3, x4]
            y = [y1, y2, y3, y4]
           
            x_new , y_new = reproject_coordinates_batch(x, y, pSpatial_reference)
            x1, x2, x3, x4 = x_new
            y1, y2, y3, y4 = y_new       
            coordinates = [(x1, y1), (x2, y2), (x3, y3), (x4, y4), (x1, y1)]


            ring = ogr.Geometry(ogr.wkbLinearRing)
            for x, y in coordinates:
                ring.AddPoint(x, y)

            pPolygon = ogr.Geometry(ogr.wkbPolygon)
            pPolygon.AddGeometry(ring)

            aCoords = np.full((5,2), -9999.0, dtype=float)
            for i, (x, y) in enumerate(coordinates):
                aCoords[i, 0] = x
                aCoords[i, 1] = y

            dummy1= np.array(aCoords)
            dLongitude_center = np.mean(aCoords[0:4,0])
            dLatitude_center = np.mean(aCoords[0:4,1])     

            iFlag = False
            if pPolygon.Within(pBoundary):
                iFlag = True
            else:
                #then check intersection
                if pPolygon.Intersects(pBoundary):
                    iFlag = True
                else:
                    pass

            if ( iFlag == True ):
                aSquare, dArea = add_cell_into_list(aSquare, lCellID, iRow, iColumn, dLongitude_center,dLatitude_center, dummy1 ) 
        
                #save feature
                pFeature.SetGeometry(pPolygon)
                pFeature.SetField("cellid", lCellID)
                pFeature.SetField("longitude", dLongitude_center )
                pFeature.SetField("latitude", dLatitude_center )
                pFeature.SetField("area", dArea )
                pLayer.CreateFeature(pFeature)

                pass

    
    aSquare_out = list()
    aSquare_middle = list()
 
    ncell = len(aSquare)
    aCellID  = list()
    for i in range(ncell):
        pCell = aSquare[i]
        lCellID = pCell.lCellID
        aCellID.append(lCellID)

    for i in range(ncell):
        pCell = aSquare[i]
        aNeighbor_land = pCell.aNeighbor_land   #including both holes and maps land cutoff by boundary
        nNeighbor_land = pCell.nNeighbor
        aNeighbor_land_update = list()
        aNeighbor_land_virtual = list()
        nNeighbor_land_update = 0 
        for j in range(nNeighbor_land): #loop all land neighbors
            lNeighbor = int(aNeighbor_land[j])
            if lNeighbor in aCellID:
                nNeighbor_land_update = nNeighbor_land_update + 1 
                aNeighbor_land_update.append(lNeighbor)
            else:
                #a hole or boundary mpas land cell
                aNeighbor_land_virtual.append(lNeighbor)
                
        pCell.nNeighbor= len(aNeighbor_land_update)
        pCell.aNeighbor = aNeighbor_land_update        
        pCell.aNeighbor_land = aNeighbor_land_update
        pCell.nNeighbor_land= len(aNeighbor_land_update)   
        pCell.aNeighbor_land_virtual = aNeighbor_land_virtual   
        
        pCell.nNeighbor_land_virtual = len(aNeighbor_land_virtual)
        aSquare_middle.append(pCell)

    #add hole back
    for i in range(ncell):
        pCell = aSquare_middle[i]  
          
        if pCell.nNeighbor_land_virtual ==1:  #only one virtual land means it is likely next to a hole 
            lNeighbor_hole = pCell.aNeighbor_land_virtual[0]
            #now find its row and column indices
            #id start with 1 so we need to refind the row and column index
            iRow, iColumn = index_to_row_col(lNeighbor_hole, ncolumn_in)          
            lCellID = (iRow-1) * ncolumn_in + iColumn
            if lCellID != lNeighbor_hole:
                print("error")
                return
            
            #now build the cell    
            #define a polygon here
            x1 = xleft + ((iColumn-1) * xspacing)
            y1 = ybottom + ((iRow-1) * yspacing)

            x2 = xleft + ((iColumn ) * xspacing)
            y2 = ybottom + ((iRow-1) * yspacing)     

            x3 = xleft + ((iColumn ) * xspacing)
            y3 = ybottom + ((iRow ) * yspacing)

            x4 = xleft + ((iColumn-1) * xspacing)
            y4 = ybottom + ((iRow ) * yspacing)  

            x = [x1, x2, x3, x4]
            y = [y1, y2, y3, y4]
           
            x_new , y_new = reproject_coordinates_batch(x, y, pSpatial_reference)
            x1, x2, x3, x4 = x_new
            y1, y2, y3, y4 = y_new       
            coordinates = [(x1, y1), (x2, y2), (x3, y3), (x4, y4), (x1, y1)]        

            ring = ogr.Geometry(ogr.wkbLinearRing)
            for x, y in coordinates:
                ring.AddPoint(x, y)

            pPolygon = ogr.Geometry(ogr.wkbPolygon)
            pPolygon.AddGeometry(ring)

            aCoords = np.full((5,2), -9999.0, dtype=float)
            for i, (x, y) in enumerate(coordinates):
                aCoords[i, 0] = x
                aCoords[i, 1] = y
    

            dummy1= np.array(aCoords)
            dLongitude_center = np.mean(aCoords[0:4,0])
            dLatitude_center = np.mean(aCoords[0:4,1])   
            
            if lCellID not in aCellID:    
                aSquare_middle, dArea = add_cell_into_list(aSquare_middle, lCellID, iRow, iColumn, dLongitude_center,dLatitude_center, dummy1 )        
                aCellID.append(lCellID)

                pFeature.SetGeometry(pPolygon)
                pFeature.SetField("cellid", int(lCellID) )
                pFeature.SetField("longitude", dLongitude_center )
                pFeature.SetField("latitude", dLatitude_center )
                pFeature.SetField("area", dArea )                
                pLayer.CreateFeature(pFeature)

            else:
                #this hole was added already, but we need to update the neighbor information
                pCell.aNeighbor_land.append(lCellID)
                pCell.nNeighbor_land = pCell.nNeighbor_land + 1
                pCell.aNeighbor_land_virtual = None
                pCell.nNeighbor_land_virtual = 0
                pass

    #update
    ncell = len(aSquare_middle)
    for i in range(ncell):
        pCell = aSquare_middle[i]
        aNeighbor_land_update = list()   
        aNeighbor_land = pCell.aNeighbor_land                    
        nNeighbor_land = pCell.nNeighbor_land
        aNeighbor_land_virtual_update = list()      
        aNeighbor_land_virtual = pCell.aNeighbor_land_virtual
        nNeighbor_land_virtual = pCell.nNeighbor_land_virtual       
        
        for j in range(nNeighbor_land):
            lNeighbor = int(aNeighbor_land[j])
            
            if lNeighbor in aCellID:
                aNeighbor_land_update.append(lNeighbor)
                
                pass
            else:
                #this is a land cell in mpas, but it may be clipped by boundary
                pass

        #for book keeping only        
        for j in range(nNeighbor_land_virtual):
            lNeighbor = int(aNeighbor_land_virtual[j])
            if lNeighbor in aCellID:
                #this cell is actually not virtual anymore                    
                aNeighbor_land_update.append(lNeighbor)
            else:
                aNeighbor_land_virtual_update.append(lNeighbor)
   
        pCell.nNeighbor= len(aNeighbor_land_update)
        pCell.aNeighbor = aNeighbor_land_update        
        pCell.aNeighbor_land = aNeighbor_land_update
        pCell.nNeighbor_land= len(aNeighbor_land_update)   
        pCell.aNeighbor_land_virtual = aNeighbor_land_virtual_update   #for book keeping only
        pCell.nNeighbor_land_virtual = len(aNeighbor_land_virtual_update)
        
        aSquare_out.append(pCell)

    for pSquare in aSquare_out:
        aNeighbor = pSquare.aNeighbor
        pSquare.aNeighbor_distance=list()
        for lCellID1 in aNeighbor:
            for pSquare1 in aSquare_out:
                if pSquare1.lCellID == lCellID1:
                    dDistance = pSquare.pVertex_center.calculate_distance( pSquare1.pVertex_center )
                    pSquare.aNeighbor_distance.append(dDistance)
                    break

    pDataset = pLayer = pFeature  = None  
    return aSquare_out

