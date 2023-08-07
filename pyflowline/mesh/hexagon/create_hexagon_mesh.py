#create a rectangle latitude/longitude based mesh
#we will use some GIS way to define it
#longitude left and latitude bottom and nrow and ncolumn and resolution is used to define the rectangle
#because it is mesh, it represent the edge instead of center
#we will use gdal api for most operations
import os
import numpy as np
from osgeo import ogr, osr
from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_cell
from pyflowline.algorithms.auxiliary.find_index_in_list import check_if_duplicates
from pyflowline.external.pyearth.gis.gdal.gdal_functions import reproject_coordinates_batch

def index_to_row_col(index, num_columns):
    index -= 1  # Adjust for 1-based indexing
    row = index // num_columns + 1
    col = index % num_columns + 1
    return row, col

def create_hexagon_mesh(iFlag_rotation_in, 
        dX_left_in, 
        dY_bot_in, 
        dResolution_meter_in, 
        ncolumn_in, 
        nrow_in, 
        sFilename_output_in, 
        sFilename_spatial_reference_in,
        pBoundary_in):
    """
    _summary_

    Args:
        iFlag_rotation_in (_type_): _description_
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
    
    if os.path.exists(sFilename_spatial_reference_in): 
        pass
    else:
        print('Spatial reference file is missing')
        return

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
    ybottom = dY_bot_in
    dArea = np.power(dResolution_meter_in,2.0)
    #hexagon edge
    dLength_edge = np.sqrt(  2.0 * dArea / (3.0* np.sqrt(3.0))  )
    dLength_half_edge = 0.5 * dLength_edge

    def add_cell_into_list1(aList, lCellID, iRow, iColumn, dLongitude_center, dLatitude_center, aCoords ):
        pHexagon = convert_gcs_coordinates_to_cell(1, dLongitude_center, dLatitude_center, aCoords)
        pHexagon.lCellID = lCellID
        dArea = pHexagon.calculate_cell_area()
        pHexagon.calculate_edge_length()
        
        lCellID_center = lCellID
        #build topoloy
        aNeighbor=list()
        if iRow > 1 : #1 and 2
            if iRow %2 ==1:
                if iColumn !=1:
                    iRow_dummy = iRow - 1
                    iColumn_dummy = iColumn - 1
                    lCellID1 = (iRow_dummy-1) * ncolumn_in + iColumn_dummy  
                    aNeighbor.append(lCellID1)
                    
                iRow_dummy = iRow - 1
                iColumn_dummy = iColumn 
                lCellID2 = (iRow_dummy-1) * ncolumn_in + iColumn_dummy
                aNeighbor.append(lCellID2)
                
            else:
                iRow_dummy = iRow - 1
                iColumn_dummy = iColumn 
                lCellID1 = (iRow_dummy-1) * ncolumn_in + iColumn_dummy
                aNeighbor.append(lCellID1)
                if iColumn!=ncolumn_in:
                    iRow_dummy = iRow - 1
                    iColumn_dummy = iColumn + 1
                    lCellID2 = (iRow_dummy-1) * ncolumn_in + iColumn_dummy
                    
                    aNeighbor.append(lCellID2)

        if iColumn < ncolumn_in:#3
            iRow_dummy = iRow 
            iColumn_dummy = iColumn + 1
            lCellID3 = (iRow_dummy-1) * ncolumn_in + iColumn_dummy
            aNeighbor.append(lCellID3)  

        
        if iRow < nrow_in :#4 and 5
            if iRow %2 ==0:
                if iColumn!=ncolumn_in:
                    iRow_dummy = iRow + 1
                    iColumn_dummy = iColumn + 1
                    lCellID4 = (iRow_dummy-1) * ncolumn_in + iColumn_dummy
                    aNeighbor.append(lCellID4)

                iRow_dummy = iRow + 1
                iColumn_dummy = iColumn 
                lCellID5 = (iRow_dummy-1) * ncolumn_in + iColumn_dummy
                aNeighbor.append(lCellID5)
                
            else:
                iRow_dummy = iRow + 1
                iColumn_dummy = iColumn 
                lCellID4 = (iRow_dummy-1) * ncolumn_in + iColumn_dummy
                aNeighbor.append(lCellID4)
                if iColumn != 1:
                    iRow_dummy = iRow + 1
                    iColumn_dummy = iColumn -1
                    lCellID5 = (iRow_dummy-1) * ncolumn_in + iColumn_dummy
                    aNeighbor.append(lCellID5)
                   
        if iColumn > 1:#6
            iRow_dummy = iRow 
            iColumn_dummy = iColumn - 1
            lCellID6 = (iRow_dummy-1) * ncolumn_in + iColumn_dummy
            aNeighbor.append(lCellID6)

        if check_if_duplicates(aNeighbor) == 0:
            print('error')  

        pHexagon.aNeighbor = aNeighbor
        pHexagon.nNeighbor = len(aNeighbor)
        pHexagon.aNeighbor_land= aNeighbor
        pHexagon.nNeighbor_land= pHexagon.nNeighbor
        aList.append(pHexagon)  
        return aList, dArea

    def add_cell_into_list2(aList, lCellID, iRow, iColumn, dLongitude_center, dLatitude_center, aCoords ):
        pHexagon = convert_gcs_coordinates_to_cell(1, dLongitude_center, dLatitude_center, aCoords)
        pHexagon.lCellID = lCellID
        dArea = pHexagon.calculate_cell_area()
        pHexagon.dArea = dArea
        pHexagon.calculate_edge_length() 
        
                      
        aNeighbor=list()
        if iRow > 1:#1 under
            #lCellID0 = lCellID_center - 1
            iRow_dummy = iRow -1
            iColumn_dummy = iColumn 
            lCellID1 = (iRow_dummy-1) * ncolumn_in + iColumn_dummy
            aNeighbor.append(lCellID1)
        
        if iColumn  < ncolumn_in  : #2 and 3
            if iColumn %2 ==1:
                iRow_dummy = iRow 
                iColumn_dummy = iColumn + 1
                lCellID2 = (iRow_dummy-1) * ncolumn_in + iColumn_dummy 
                aNeighbor.append(lCellID2)
                if iRow !=nrow_in:
                    iRow_dummy = iRow + 1
                    iColumn_dummy = iColumn + 1
                    lCellID3 = (iRow_dummy-1) * ncolumn_in + iColumn_dummy 
                    aNeighbor.append(lCellID3)
            else:
                if iRow!=1:
                    iRow_dummy = iRow 
                    iColumn_dummy = iColumn + 1
                    lCellID2 = (iRow_dummy-1) * ncolumn_in + iColumn_dummy 
                    aNeighbor.append(lCellID2)

                iRow_dummy = iRow + 1
                iColumn_dummy = iColumn + 1
                lCellID3 = (iRow_dummy-1) * ncolumn_in + iColumn_dummy 
                aNeighbor.append(lCellID3)
                

        if iRow < nrow_in:#4
            iRow_dummy = iRow + 1
            iColumn_dummy = iColumn 
            lCellID4 = (iRow_dummy-1) * ncolumn_in + iColumn_dummy
            aNeighbor.append(lCellID4)

        if iColumn> 1:#5 and 6
            if iColumn %2 ==1:
                if iRow != nrow_in:
                    iRow_dummy = iRow + 1
                    iColumn_dummy = iColumn -1
                    lCellID5 = (iRow_dummy-1) * ncolumn_in + iColumn_dummy
                    aNeighbor.append(lCellID5)

                iRow_dummy = iRow 
                iColumn_dummy = iColumn -1
                lCellID6 = (iRow_dummy-1) * ncolumn_in + iColumn_dummy
                aNeighbor.append(lCellID6)
            else:
                iRow_dummy = iRow 
                iColumn_dummy = iColumn - 1
                lCellID5 = (iRow_dummy-1) * ncolumn_in + iColumn_dummy
                aNeighbor.append(lCellID5)
                if iRow!=1:
                    iRow_dummy = iRow - 1 
                    iColumn_dummy = iColumn - 1
                    lCellID6 = (iRow_dummy-1) * ncolumn_in + iColumn_dummy
                    aNeighbor.append(lCellID6)
        
        
        if check_if_duplicates(aNeighbor) == 0:
            print('error')  

        pHexagon.aNeighbor = aNeighbor
        pHexagon.nNeighbor = len(aNeighbor)
        pHexagon.aNeighbor_land= aNeighbor
        pHexagon.nNeighbor_land= pHexagon.nNeighbor
        aList.append(pHexagon)
        return aList, dArea

    #geojson
    aHexagon=list()
    #.........
    #change the order because mpas uses counter-clock wise to store the vertices
    #we will also start from the lower-left corner, and then go to the right and then go up
    #so the final index will be like this
    #3 4
    #1 2    
    #iFlag_rotation_in = 0, easy to row-column order
    #---------(x4,y4)
    #--(x5,y5)--------(x3,y3)
    #--(x6,y6)--------(x2,y2)
    #---------(x1,y1)
    #iFlag_rotation_in = 1, easy to column-row order
    #--------(x5,y6)----(x4,y4)
    #--(x6,y6)----------------(x3,y3)
    #--------(x1,y1)----(x2,y2)
    #...............#
    if iFlag_rotation_in == 0:

        dX_shift = dLength_half_edge * np.sqrt(3.0)
        dY_shift = dLength_half_edge
        dX_spacing = dLength_edge * np.sqrt(3.0)
        dY_spacing = dLength_edge * 1.5
        for iRow in range(1, nrow_in+1):
            for iColumn in range(1, ncolumn_in+1):
                #using global id to identify the cell
                lCellID = (iRow-1) * ncolumn_in + iColumn
                if iRow % 2 == 1 : #odd
                #define a polygon here
                    x1 = xleft + (iColumn-1) * dX_spacing
                    y1 = ybottom + (iRow -1) * dY_spacing
                else:
                    x1 = xleft + (iColumn-1) * dX_spacing + dX_shift
                    y1 = ybottom + (iRow -1) * dY_spacing    

                x2 = x1 + dX_shift
                y2 = y1 + dY_shift
    
                x3 = x2
                y3 = y2 + dLength_edge            
    
                x4 = x1 
                y4 = y2 + dY_spacing
    
                x5 = x1- dX_shift 
                y5 = y3 
    
                x6 = x5
                y6 = y2                  
                
                x = [x1, x2, x3, x4, x5, x6]
                y = [y1, y2, y3, y4, y5, y6]

                x_new , y_new = reproject_coordinates_batch(x, y, pSpatial_reference, \
                    spatial_reference_target = pSpatial_reference_gcs)

                x1, x2, x3, x4, x5, x6 = x_new
                y1, y2, y3, y4, y5, y6 = y_new
                coordinates = [(x1, y1), (x2, y2), (x3, y3), (x4, y4), (x5, y5), (x6, y6), (x1, y1)]
    
                ring = ogr.Geometry(ogr.wkbLinearRing)
                for x, y in coordinates:
                    ring.AddPoint(x, y)
                
                pPolygon = ogr.Geometry(ogr.wkbPolygon)
                pPolygon.AddGeometry(ring)    

                aCoords = np.full((7,2), -9999.0, dtype=float)
                for i, (x, y) in enumerate(coordinates):
                    aCoords[i, 0] = x
                    aCoords[i, 1] = y
           
                dummy1= np.array(aCoords)
                dLongitude_center = np.mean(aCoords[0:6,0])
                dLatitude_center = np.mean(aCoords[0:6,1])                

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
                    aHexagon, dArea = add_cell_into_list1(aHexagon, lCellID, iRow, iColumn, dLongitude_center,dLatitude_center, dummy1 ) 
        
                    #save feature
                    pFeature.SetGeometry(pPolygon)
                    pFeature.SetField("cellid", lCellID) 
                    pFeature.SetField("longitude", dLongitude_center )
                    pFeature.SetField("latitude", dLatitude_center )
                    pFeature.SetField("area", dArea )
                    pLayer.CreateFeature(pFeature)
  
                    pass
                else:
                    #this cell center is out of boundary
                    continue
              
        pass
    else:
        dX_shift = 0.5 * dLength_edge #short displacement 
        dY_shift = 0.5 * dLength_edge * np.sqrt(3.0) #long displacement
        dX_spacing = dLength_edge * 1.5 #cell to cell displacement in x
        dY_spacing = dLength_edge * np.sqrt(3.0) #cell to cell displacement in y
        for iColumn in range(1, ncolumn_in+1):
            for iRow in range(1, nrow_in+1):
                if iColumn % 2 == 0 :
                #define a polygon here
                    x1 = xleft + (iColumn-1) * dX_spacing
                    y1 = ybottom + (iRow-2) * dY_spacing + dY_shift
                else:
                    x1 = xleft + (iColumn-1) * dX_spacing 
                    y1 = ybottom + (iRow-1) * dY_spacing     
    
                x2 = x1 + dLength_edge
                y2 = y1  

                x3 = x2 + dX_shift
                y3 = y2 + dY_shift
    
                x4 = x1 + dLength_edge
                y4 = y1 + dY_shift * 2.0

                x5 = x1 
                y5 = y1 + dY_shift * 2.0
    
                x6 = x1 - dX_shift
                y6 = y1 + dY_shift                          
               
                x = [x1, x2, x3, x4, x5, x6]
                y = [y1, y2, y3, y4, y5, y6]

                x_new , y_new = reproject_coordinates_batch(x, y, pSpatial_reference)
                x1, x2, x3, x4, x5, x6 = x_new
                y1, y2, y3, y4, y5, y6 = y_new
                coordinates = [(x1, y1), (x2, y2), (x3, y3), (x4, y4), (x5, y5), (x6, y6), (x1, y1)]

                ring = ogr.Geometry(ogr.wkbLinearRing)
                
                for x, y in coordinates:
                    ring.AddPoint(x, y)

                pPolygon = ogr.Geometry(ogr.wkbPolygon)
                pPolygon.AddGeometry(ring)

                aCoords = np.full((7,2), -9999.0, dtype=float)
                for i, (x, y) in enumerate(coordinates):
                    aCoords[i, 0] = x
                    aCoords[i, 1] = y
                    
                dummy1= np.array(aCoords)
                dLongitude_center = np.mean(aCoords[0:6,0])
                dLatitude_center = np.mean(aCoords[0:6,1])     
                iFlag == False
                if pPolygon.Within(pBoundary):
                    iFlag = True
                else:
                    #then check intersection
                    if pPolygon.Intersects(pBoundary):
                        iFlag = True
                    else:
                        pass
                if ( iFlag == True ):  
                    aHexagon, dArea = add_cell_into_list2(aHexagon, lCellID, iRow, iColumn, 
                                                          dLongitude_center,dLatitude_center, 
                                                          dummy1 )
                    
                    pFeature.SetGeometry(pPolygon)
                    pFeature.SetField("cellid", lCellID)
                    pFeature.SetField("longitude", dLongitude_center )
                    pFeature.SetField("latitude", dLatitude_center )
                    pFeature.SetField("area", dArea )
                    pLayer.CreateFeature(pFeature)

                    pass
                else:
                    #out of bound
                    pass
       
        
         

    #maybe rebuild topology?
    
    aHexagon_out = list()
    aHexagon_middle = list()
    ncell = len(aHexagon)
    aCellID  = list()
    for i in range(ncell):
        pCell = aHexagon[i]
        lCellID = pCell.lCellID
        aCellID.append(lCellID)

    for i in range(ncell):
        pCell = aHexagon[i]
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
        aHexagon_middle.append(pCell)

    #add hole back
    for i in range(ncell):
        pCell = aHexagon_middle[i]  
          
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
            if iFlag_rotation_in == 0:
                dX_shift = dLength_half_edge * np.sqrt(3.0)
                dY_shift = dLength_half_edge
                dX_spacing = dLength_edge * np.sqrt(3.0)
                dY_spacing = dLength_edge * 1.5
                if iRow % 2 == 1 : #odd
                #define a polygon here
                    x1 = xleft + (iColumn-1) * dX_spacing
                    y1 = ybottom + (iRow -1) * dY_spacing
                else:
                    x1 = xleft + (iColumn-1) * dX_spacing + dX_shift
                    y1 = ybottom + (iRow -1) * dY_spacing    

                x2 = x1 + dX_shift
                y2 = y1 + dY_shift
    
                x3 = x2
                y3 = y2 + dLength_edge            
    
                x4 = x1 
                y4 = y2 + dY_spacing
    
                x5 = x1- dX_shift 
                y5 = y3 
    
                x6 = x5
                y6 = y2                  
                
                x = [x1, x2, x3, x4, x5, x6]
                y = [y1, y2, y3, y4, y5, y6]

                x_new , y_new = reproject_coordinates_batch(x, y, pSpatial_reference, \
                    spatial_reference_target = pSpatial_reference_gcs)

                x1, x2, x3, x4, x5, x6 = x_new
                y1, y2, y3, y4, y5, y6 = y_new
                coordinates = [(x1, y1), (x2, y2), (x3, y3), (x4, y4), (x5, y5), (x6, y6), (x1, y1)]
    
                ring = ogr.Geometry(ogr.wkbLinearRing)
                for x, y in coordinates:
                    ring.AddPoint(x, y)
                
                pPolygon = ogr.Geometry(ogr.wkbPolygon)
                pPolygon.AddGeometry(ring)    

                aCoords = np.full((7,2), -9999.0, dtype=float)
                for i, (x, y) in enumerate(coordinates):
                    aCoords[i, 0] = x
                    aCoords[i, 1] = y
           
                dummy1= np.array(aCoords)
                dLongitude_center = np.mean(aCoords[0:6,0])
                dLatitude_center = np.mean(aCoords[0:6,1])   
                

                dummy1= np.array(aCoords)
                dLongitude_center = np.mean(aCoords[0:4,0])
                dLatitude_center = np.mean(aCoords[0:4,1])   

                if lCellID not in aCellID:    
                    aHexagon_middle, dArea = add_cell_into_list1(aHexagon_middle, lCellID, iRow, iColumn, dLongitude_center,dLatitude_center, dummy1 )        
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
            
            else:
                dX_shift = 0.5 * dLength_edge #short displacement 
                dY_shift = 0.5 * dLength_edge * np.sqrt(3.0) #long displacement
                dX_spacing = dLength_edge * 1.5 #cell to cell displacement in x
                dY_spacing = dLength_edge * np.sqrt(3.0) #cell to cell displacement in y
                if iColumn % 2 == 0 :
                #define a polygon here
                    x1 = xleft + (iColumn-1) * dX_spacing
                    y1 = ybottom + (iRow-2) * dY_spacing + dY_shift
                else:
                    x1 = xleft + (iColumn-1) * dX_spacing 
                    y1 = ybottom + (iRow-1) * dY_spacing     
    
                x2 = x1 + dLength_edge
                y2 = y1  

                x3 = x2 + dX_shift
                y3 = y2 + dY_shift
    
                x4 = x1 + dLength_edge
                y4 = y1 + dY_shift * 2.0

                x5 = x1 
                y5 = y1 + dY_shift * 2.0
    
                x6 = x1 - dX_shift
                y6 = y1 + dY_shift      

                x_new , y_new = reproject_coordinates_batch(x, y, pSpatial_reference, \
                    spatial_reference_target = pSpatial_reference_gcs)

                x1, x2, x3, x4, x5, x6 = x_new
                y1, y2, y3, y4, y5, y6 = y_new
                coordinates = [(x1, y1), (x2, y2), (x3, y3), (x4, y4), (x5, y5), (x6, y6), (x1, y1)]
    
                ring = ogr.Geometry(ogr.wkbLinearRing)
                for x, y in coordinates:
                    ring.AddPoint(x, y)
                
                pPolygon = ogr.Geometry(ogr.wkbPolygon)
                pPolygon.AddGeometry(ring)    

                aCoords = np.full((7,2), -9999.0, dtype=float)
                for i, (x, y) in enumerate(coordinates):
                    aCoords[i, 0] = x
                    aCoords[i, 1] = y
           
                dummy1= np.array(aCoords)
                dLongitude_center = np.mean(aCoords[0:6,0])
                dLatitude_center = np.mean(aCoords[0:6,1])   
                

                dummy1= np.array(aCoords)
                dLongitude_center = np.mean(aCoords[0:4,0])
                dLatitude_center = np.mean(aCoords[0:4,1])   

                if lCellID not in aCellID:    
                    aHexagon_middle, dArea = add_cell_into_list1(aHexagon_middle, lCellID, iRow, iColumn, dLongitude_center,dLatitude_center, dummy1 )        
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
                pass

    ncell = len(aHexagon_middle)
    for i in range(ncell):
        pCell = aHexagon_middle[i]
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
        
        aHexagon_out.append(pCell)
        
    #calculate neighbor distance
    for pHexagon in aHexagon_out:
        aNeighbor = pHexagon.aNeighbor
        pHexagon.aNeighbor_distance=list()
        for lCellID1 in aNeighbor:
            for pHexagon1 in aHexagon_out:
                if pHexagon1.lCellID == lCellID1:
                    dDistance = pHexagon.pVertex_center.calculate_distance( pHexagon1.pVertex_center )
                    pHexagon.aNeighbor_distance.append(dDistance)
                    break

    pDataset = pLayer = pFeature  = None 
    return aHexagon_out



