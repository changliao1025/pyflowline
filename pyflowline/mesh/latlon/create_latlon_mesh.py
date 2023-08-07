#create a rectangle latitude/longitude based mesh
#we will use some GIS way to define it
#longitude left and latitude bottom and nrow and ncolumn and resolution is used to define the rectangle
#because it is mesh, it represent the edge instead of center
#we will use gdal api for most operations
import os
import numpy as np
from osgeo import ogr, osr
from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_cell
def create_latlon_mesh(dLongitude_left_in, 
                       dLatitude_bot_in, 
                       dResolution_degree_in, 
                       ncolumn_in, nrow_in, 
                       sFilename_output_in,
                       pBoundary_in):
    """
    _summary_

    Args:
        dLongitude_left_in (_type_): _description_
        dLatitude_bot_in (_type_): _description_
        dResolution_degree_in (_type_): _description_
        ncolumn_in (_type_): _description_
        nrow_in (_type_): _description_
        sFilename_output_in (_type_): _description_

    Returns:
        _type_: _description_
    """
    #for the reason that a geometry object will be crash if the associated dataset is closed, we must pass wkt string
    #https://gdal.org/api/python_gotchas.html
    pBoundary = ogr.CreateGeometryFromWkt(pBoundary_in)   
    
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
    pLayer.CreateField(ogr.FieldDefn('cellid', ogr.OFTInteger64)) #long type for high resolution
    pLayer.CreateField(ogr.FieldDefn('longitude', ogr.OFTReal)) #long type for high resolution
    pLayer.CreateField(ogr.FieldDefn('latitude', ogr.OFTReal)) #long type for high resolution
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
    aLatlon = list()
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

            coordinates = [(x1, y1), (x2, y2), (x3, y3), (x4, y4), (x1, y1)]

            ring = ogr.Geometry(ogr.wkbLinearRing)
            for x, y in coordinates:
                ring.AddPoint(x, y)

            pPolygon = ogr.Geometry(ogr.wkbPolygon)
            pPolygon.AddGeometry(ring)

            dLongitude_center = (x1 + x2 + x3 + x4)/4.0
            dLatitude_center = (y1 + y2 + y3 + y4)/4.0
            aCoords = np.full((5,2), -9999.0, dtype=float)
            for i, (x, y) in enumerate(coordinates):
                aCoords[i, 0] = x
                aCoords[i, 1] = y         

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
                dummy1= np.array(aCoords)           
                pLatlon = convert_gcs_coordinates_to_cell(3, dLongitude_center, dLatitude_center, dummy1)
                pLatlon.lCellID = lCellID
                dArea = pLatlon.calculate_cell_area()
                pLatlon.calculate_edge_length()      
                #lCellID_center = lCellID

                aNeighbor = list()

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
                    

                pLatlon.aNeighbor = aNeighbor
                pLatlon.nNeighbor = len(aNeighbor)
                pLatlon.aNeighbor_land= aNeighbor
                pLatlon.nNeighbor_land= pLatlon.nNeighbor            
                aLatlon.append(pLatlon)

                #save feature
                pFeature.SetGeometry(pPolygon)
                pFeature.SetField("cellid", lCellID)
                pFeature.SetField("longitude", dLongitude_center )
                pFeature.SetField("latitude", dLatitude_center )
                pFeature.SetField("area", dArea )
                pLayer.CreateFeature(pFeature)

                

                pass
        
    pDataset = pLayer = pFeature  = None      

    aLatlon_out = list()
    ncell = len(aLatlon)
    aCellID  = list()
    for i in range(ncell):
        pCell = aLatlon[i]
        lCellID = pCell.lCellID
        aCellID.append(lCellID)
        
    for i in range(ncell):
        pCell = aLatlon[i]
        aNeighbor = pCell.aNeighbor
        nNeighbor = pCell.nNeighbor
        aNeighbor_new = list()
        nNeighbor_new = 0 
        for j in range(nNeighbor):
            lNeighbor = int(aNeighbor[j])
            if lNeighbor in aCellID:
                nNeighbor_new = nNeighbor_new + 1 
                aNeighbor_new.append(lNeighbor)

        #for latlon, there is no ocean concept
        pCell.nNeighbor= len(aNeighbor_new)
        pCell.aNeighbor = aNeighbor_new
        pCell.nNeighbor_land= len(aNeighbor_new)
        pCell.aNeighbor_land = aNeighbor_new
        pCell.nNeighbor_ocean = pCell.nVertex - pCell.nNeighbor_land
        aLatlon_out.append(pCell)

    #calculate neighbor distance
    for pLatlon in aLatlon_out:
        aNeighbor = pLatlon.aNeighbor
        pLatlon.aNeighbor_distance=list()
        for lCellID1 in aNeighbor:
            for pLatlon1 in aLatlon_out:
                if pLatlon1.lCellID == lCellID1:
                    dDistance = pLatlon.pVertex_center.calculate_distance( pLatlon1.pVertex_center )
                    pLatlon.aNeighbor_distance.append(dDistance)
                    break

    return aLatlon_out





