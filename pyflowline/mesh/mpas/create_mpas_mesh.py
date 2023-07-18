import os
import math
import importlib
import numpy as np
from osgeo import ogr, osr, gdal
import netCDF4 as nc
from pyflowline.formats.convert_attributes import convert_gcs_attributes_to_cell
gdal.UseExceptions()  
iFlag_cython = importlib.util.find_spec("cython") 
if iFlag_cython is not None:
    from pyflowline.algorithms.cython.kernel import convert_360_to_180
else:
    from pyflowline.external.pyearth.gis.gdal.gdal_functions import convert_360_to_180

def create_mpas_mesh(iFlag_global_in, 
    iFlag_use_mesh_dem, 
    iFlag_save_mesh_in,     
    sFilename_mesh_netcdf_in, 
    sFilename_output_in,
    iFlag_antarctic_in=None,
    pBoundary_in = None): 
    """
    Create a MPAS mesh

    Args:
        iFlag_global_in (int): _description_
        iFlag_use_mesh_dem (int): _description_
        iFlag_save_mesh_in (int): _description_
        pBoundary_in (_type_): _description_
        sFilename_mesh_netcdf_in (_type_): _description_
        sFilename_output_in (_type_): _description_

    Returns:
        _type_: _description_
    """

    if iFlag_antarctic_in is None:
        iFlag_antarctic=0
    else:
        iFlag_antarctic=iFlag_antarctic_in

    if pBoundary_in is None:
        pBoundary = None
    else:
        #for the reason that a geometry object will be crash if the associated dataset is closed, we must pass wkt string
        #https://gdal.org/api/python_gotchas.html
        pBoundary = ogr.CreateGeometryFromWkt(pBoundary_in)
        #pBoundary = pBoundary_in #only works for shapely geometry object
       
    if (os.path.exists(sFilename_mesh_netcdf_in)):
        pass
    else:
        print('Mesh file does not exist!')
        return
    
    if os.path.exists(sFilename_output_in):  
        os.remove(sFilename_output_in)

    

    pDatasets_in = nc.Dataset(sFilename_mesh_netcdf_in)

    netcdf_format = pDatasets_in.file_format
    pDriver_geojson = ogr.GetDriverByName('GeoJSON')     
    pSpatial_reference_gcs = osr.SpatialReference()  
    pSpatial_reference_gcs.ImportFromEPSG(4326)    # WGS84 lat/lon    
    #geojson
    if iFlag_save_mesh_in ==1:
        pDataset = pDriver_geojson.CreateDataSource(sFilename_output_in)
        pLayer = pDataset.CreateLayer('cell', pSpatial_reference_gcs, ogr.wkbPolygon)
        # Add one attribute
        pLayer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger64)) #long type for high resolution
        pLayer.CreateField(ogr.FieldDefn('lon', ogr.OFTReal)) #long type for high resolution
        pLayer.CreateField(ogr.FieldDefn('lat', ogr.OFTReal)) #long type for high resolution
        pArea_field = ogr.FieldDefn('area', ogr.OFTReal)
        pArea_field.SetWidth(20)
        pArea_field.SetPrecision(2)
        pLayer.CreateField(pArea_field)
        if iFlag_use_mesh_dem == 1:
            pLayer.CreateField(ogr.FieldDefn('elev', ogr.OFTReal)) #float type for high resolution
            pLayer.CreateField(ogr.FieldDefn('elev0', ogr.OFTReal)) #float type for high resolution
        else:
            pass
    
        pLayerDefn = pLayer.GetLayerDefn()
        pFeature = ogr.Feature(pLayerDefn)        

    #read new netcdf
    for sKey, aValue in pDatasets_in.variables.items():      
        #we need to filter out unused grids based on mpas specs
        if sKey == 'latCell':
            latCell0 = aValue
        else:
            pass
        if sKey == 'lonCell':
            lonCell0 = aValue 
        else:
            pass
        
        if sKey == 'edgesOnCell':
            edgesOnCell0 = aValue 
        else:
            pass
            
        if sKey == 'cellsOnCell':
            cellsOnCell0 = aValue 
        else:
            pass

        if sKey == 'cellsOnEdge':
            cellsOnEdge0 = aValue 
        else:
            pass
            
        if sKey == 'verticesOnCell':
            verticesOnCell0 = aValue
        else:
            pass
        
        if sKey == 'verticesOnEdge':
            verticesOnEdge0 = aValue 
        else:
            pass

        if sKey == 'indexToCellID':
            indexToCellID0 = aValue 
        else:
            pass

        if sKey == 'indexToEdgeID':
            indexToEdgeID0 = aValue 
        else:
            pass

        if sKey == 'indexToVertexID':
            indexToVertexID0 = aValue 
        else:
            pass

        if sKey == 'lonVertex':
            lonVertex0 = aValue 
        else:
            pass

        if sKey == 'latVertex':
            latVertex0 = aValue 
        else:
            pass

        if sKey == 'areaCell':
            areaCell0 = aValue 
        else:
            pass

        if sKey == 'bed_elevation':
            bed_elevation0 = aValue 
        else:
            pass

        if sKey == 'ice_thickness':
            ice_thickness0 = aValue 
        else:
            pass

        if sKey == 'areaCell':
            areaCell0 = aValue 
        else:
            pass

        if sKey == 'dcEdge':
            dcEdge0 = aValue 
        else:
            pass

        if sKey == 'bed_elevation_profile':
            bed_elevation_profile0 = aValue 
        else:
            pass
        

    aLatitudeVertex = latVertex0[:] / math.pi * 180
    aLongitudeVertex = lonVertex0[:] / math.pi * 180    
    #convert unit 
    aLatitudeCell = latCell0[:] / math.pi * 180
    aLongitudeCell = lonCell0[:] / math.pi * 180
    aCellsOnCell = cellsOnCell0[:]
    aCellOnEdge = cellsOnEdge0[:]
    aEdgesOnCell = edgesOnCell0[:]
    aVertexOnCell = verticesOnCell0[:]
    aVertexOnEdge0 = verticesOnEdge0[:]    
    aIndexToCellID = indexToCellID0[:]
    aIndexToEdgeID = indexToEdgeID0[:]
    aIndexToVertexID = indexToVertexID0[:]
    aBed_elevation = bed_elevation0[:]
    aIce_thickness = ice_thickness0[:]
    aCellArea = areaCell0[:]
    aDcEdge = dcEdge0[:]
    aBed_elevation_profile = bed_elevation_profile0[:]  #elevation    
    ncell = len(aIndexToCellID) 
    aMpas = list()

    #add a mpas cell into a list
    def add_cell_into_list(aList, i, lCellID, dArea,dElevation_mean,dElevation_profile0, aCoords  ):
        dLat = convert_360_to_180 (aLatitudeCell[i])
        dLon = convert_360_to_180 (aLongitudeCell[i])
        #vertex
        aCellOnCellIndex = np.array(aCellsOnCell[i,:])
        aEdgesOnCellIndex = np.array(aEdgesOnCell[i,:])
        aVertexOnCellIndex = np.array(aVertexOnCell[i,:])
        dummy0 = np.where(aVertexOnCellIndex > 0)
        aVertexIndex = aVertexOnCellIndex[dummy0]
        dummy1 = np.where(aEdgesOnCellIndex > 0)
        aEdgeIndex= aEdgesOnCellIndex[dummy1]
        dummy2 = np.where(aCellOnCellIndex > 0)
        aNeighborIndex= (aCellOnCellIndex[dummy2]).astype(int)
        aVertexIndexOnEdge = np.array(aVertexOnEdge0[aEdgeIndex-1,:]).astype((int))
    
        pmpas = convert_gcs_attributes_to_cell(4, dLon, dLat, aCoords, aVertexIndex, aEdgeIndex, aVertexIndexOnEdge)               
        pmpas.dArea = dArea
        pmpas.calculate_edge_length()
        pmpas.dLength_flowline = pmpas.dLength_edge #Default
        pmpas.lCellID = lCellID
        pmpas.dElevation_mean  = dElevation_mean
        pmpas.dElevation_profile0 = dElevation_profile0
        #now setup the neighbor information
        pmpas.aNeighbor=aNeighborIndex
        pmpas.nNeighbor=len(aNeighborIndex)
        if pmpas.nNeighbor != pmpas.nVertex:
            #this cell is next to the ocean boundary
            pmpas.nNeighbor_land = pmpas.nNeighbor 
            pmpas.nNeighbor_ocean = pmpas.nVertex - pmpas.nNeighbor
            pmpas.aNeighbor_land=aNeighborIndex
            pmpas.nNeighbor_land=len(aNeighborIndex)
            print(lCellID)
            print('Warning: nNeighbor != nVertex')
        else:
            pmpas.nNeighbor_land = pmpas.nNeighbor 
            pmpas.nNeighbor_ocean = 0
            pmpas.aNeighbor_land = aNeighborIndex        
            pmpas.nNeighbor_land=len(aNeighborIndex)

        aDistance=list()
        for i in range(pmpas.nNeighbor):
            #find shared edge
            lEdgeID= aEdgeIndex[i]              
            lIndex = lEdgeID-1
            dDistance = aDcEdge[lIndex]
            aDistance.append(dDistance)
            pass

        pmpas.aNeighbor_distance = aDistance
        aList.append(pmpas)
        return aList

    if iFlag_antarctic == 1:
        iFlag_remove_ice = 0
        #if it is antarctic, we dont need the boundary
        for i in range(ncell):
            #center
            dLat = convert_360_to_180 (aLatitudeCell[i])
            dLon = convert_360_to_180 (aLongitudeCell[i])
            aVertexOnCellIndex = np.array(aVertexOnCell[i,:])
            dummy0 = np.where(aVertexOnCellIndex > 0)
            aVertexIndex = aVertexOnCellIndex[dummy0]
            aLonVertex = aLongitudeVertex[aVertexIndex-1]
            aLatVertex = aLatitudeVertex[aVertexIndex-1]
            nVertex = len(aLonVertex)
            #first check if it is within the boundary
            iFlag = False
            ring = ogr.Geometry(ogr.wkbLinearRing)
            aCoords = np.full((nVertex,2), -9999.0, dtype=float)
            for j in range(nVertex):
                x1 = convert_360_to_180(aLonVertex[j])
                y1 = aLatVertex[j] 
                ring.AddPoint(x1, y1)
                aCoords[j,0] = x1
                aCoords[j,1] = y1
                pass
            if dLat < -60:
                iFlag = True
            else:  
                iFlag = False
                pass

            if ( iFlag == True ):
                lCellID = int(aIndexToCellID[i])
                dElevation_mean = float(aBed_elevation[i])
                dElevation_profile0 = float(aBed_elevation_profile[i,0])
                dThickness_ice = float( aIce_thickness[i] )
                dArea = float(aCellArea[i])

                #then check if it is ice free
                if iFlag_remove_ice == 1:
                    if dThickness_ice > 0 :
                        continue
                    else:
                        pass
                else:
                    pass      

                #call fuction to add the cell
                aMpas = add_cell_into_list(aMpas, i, lCellID, dArea, dElevation_mean, dElevation_profile0, aCoords )
                
                #save mesh cell
                if iFlag_save_mesh_in ==1:                
                    pFeature.SetGeometry(pPolygon)
                    pFeature.SetField("id", int(lCellID) )
                    pFeature.SetField("lon", dLon )
                    pFeature.SetField("lat", dLat )
                    pFeature.SetField("area", dArea )
                    if iFlag_use_mesh_dem == 1:
                        pFeature.SetField("elev", dElevation_mean )
                        pFeature.SetField("elev0", dElevation_profile0 )

                    pLayer.CreateFeature(pFeature)
            
            

    else:
        iFlag_remove_ice = 1
        for i in range(ncell):
            #center
            dLat = convert_360_to_180 (aLatitudeCell[i])
            dLon = convert_360_to_180 (aLongitudeCell[i])
            #vertex
            aVertexOnCellIndex = np.array(aVertexOnCell[i,:])
            dummy0 = np.where(aVertexOnCellIndex > 0)
            aVertexIndex = aVertexOnCellIndex[dummy0]
            aLonVertex = aLongitudeVertex[aVertexIndex-1]
            aLatVertex = aLatitudeVertex[aVertexIndex-1]
            nVertex = len(aLonVertex)
            #first check if it is within the boundary
            iFlag = False
            ring = ogr.Geometry(ogr.wkbLinearRing)
            aCoords = np.full((nVertex,2), -9999.0, dtype=float)
            for j in range(nVertex):
                x1 = convert_360_to_180(aLonVertex[j])
                y1 = aLatVertex[j] 
                ring.AddPoint(x1, y1)
                aCoords[j,0] = x1
                aCoords[j,1] = y1
                pass
            
            x1 = convert_360_to_180(aLonVertex[0])
            y1 = aLatVertex[0]
            ring.AddPoint(x1, y1) #double check            
            pPolygon = ogr.Geometry(ogr.wkbPolygon)
            pPolygon.AddGeometry(ring)
            #check within first
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
                lCellID = int(aIndexToCellID[i])
                dElevation_mean = float(aBed_elevation[i])
                dElevation_profile0 = float(aBed_elevation_profile[i,0])
                dThickness_ice = float( aIce_thickness[i] )
                dArea = float(aCellArea[i])

                #then check if it is ice free
                if iFlag_remove_ice == 1:
                    if dThickness_ice > 0 :
                        continue
                    else:
                        pass
                else:
                    pass           
                #call fuction to add the cell
                aMpas = add_cell_into_list(aMpas, i, lCellID, dArea, dElevation_mean, dElevation_profile0, aCoords )
                
                #save mesh cell
                if iFlag_save_mesh_in ==1:                
                    pFeature.SetGeometry(pPolygon)
                    pFeature.SetField("id", int(lCellID) )
                    pFeature.SetField("lon", dLon )
                    pFeature.SetField("lat", dLat )
                    pFeature.SetField("area", dArea )
                    if iFlag_use_mesh_dem == 1:
                        pFeature.SetField("elev", dElevation_mean )
                        pFeature.SetField("elev0", dElevation_profile0 )

                    pLayer.CreateFeature(pFeature)


    #for maps we need to clean some cell because they were not actually in the domain
    #besides, we need to add some smal holes back
    #to do this, we need two steps.

    if iFlag_global_in == 1:
        aMpas_out = aMpas
    else:
        aMpas_middle = list()
        aMpas_out = list()
        ncell = len(aMpas)
        #generate the list of cell ID
        aCellID  = list()
        for i in range(ncell):
            pCell = aMpas[i]
            lCellID = pCell.lCellID
            aCellID.append(lCellID)

        #first update neighbor information because some cell should have vitual land neighbor (not present in the mesh)
        for i in range(ncell):
            pCell = aMpas[i]
            aNeighbor_land = pCell.aNeighbor_land
            nNeighbor_land = pCell.nNeighbor
            aNeighbor_land_update = list()
            aNeighbor_land_virtual = list()
            nNeighbor_land_update = 0 
            for j in range(nNeighbor_land):
                lNeighbor = int(aNeighbor_land[j])
                if lNeighbor in aCellID:
                    nNeighbor_land_update = nNeighbor_land_update + 1 
                    aNeighbor_land_update.append(lNeighbor)
                else:
                    aNeighbor_land_virtual.append(lNeighbor)
                    

            pCell.aNeighbor_land = aNeighbor_land_update
            pCell.nNeighbor_land= len(aNeighbor_land_update)   
            pCell.aNeighbor_land_virtual = aNeighbor_land_virtual   
            pCell.nNeighbor_land_virtual = len(aNeighbor_land_virtual)
            aMpas_middle.append(pCell)

        #now add back small holes   
        for i in range(ncell):
            pCell = aMpas_middle[i]  
              
            if pCell.nNeighbor_land_virtual ==1:   
                lNeighbor_hole = pCell.aNeighbor_land_virtual[0]
                j = lNeighbor_hole-1
                dLat = convert_360_to_180 (aLatitudeCell[j])
                dLon = convert_360_to_180 (aLongitudeCell[j])

                #vertex
                aVertexOnCellIndex = np.array(aVertexOnCell[j,:])
                dummy0 = np.where(aVertexOnCellIndex > 0)
                aVertexIndex = aVertexOnCellIndex[dummy0]
                aLonVertex = aLongitudeVertex[aVertexIndex-1]
                aLatVertex = aLatitudeVertex[aVertexIndex-1]
                nVertex = len(aLonVertex)

                #first check if it is within the boundary
    
                ring = ogr.Geometry(ogr.wkbLinearRing)
                aCoords = np.full((nVertex,2), -9999.0, dtype=float)
                for k in range(nVertex):
                    x1 = convert_360_to_180(aLonVertex[k])
                    y1 = aLatVertex[k] 
                    ring.AddPoint(x1, y1)
                    aCoords[k,0] = x1
                    aCoords[k,1] = y1
                    pass

                lCellID = int(aIndexToCellID[j])
                dElevation_mean = float(aBed_elevation[j])
                dElevation_profile0 = float(aBed_elevation_profile[j,0])
                dArea = float(aCellArea[j])

                if lCellID not in aCellID:
                    aMpas_middle = add_cell_into_list(aMpas_middle, j, lCellID, dArea, dElevation_mean, dElevation_profile0, aCoords )
                    aCellID.append(lCellID)

                    #now we need to update the neightboring information as well
                    pCell.aNeighbor_land.append(lCellID)
                    pCell.nNeighbor_land = pCell.nNeighbor_land + 1

                    pCell.aNeighbor_land_virtual = None
                    pCell.nNeighbor_land_virtual = 0

                    #save mesh cell
                    x1 = convert_360_to_180(aLonVertex[0])
                    y1 = aLatVertex[0]
                    ring.AddPoint(x1, y1) #double check            
                    pPolygon = ogr.Geometry(ogr.wkbPolygon)
                    pPolygon.AddGeometry(ring)
                    if iFlag_save_mesh_in ==1:                
                        pFeature.SetGeometry(pPolygon)
                        pFeature.SetField("id", int(lCellID) )
                        pFeature.SetField("lon", dLon )
                        pFeature.SetField("lat", dLat )
                        pFeature.SetField("area", dArea )
                        if iFlag_use_mesh_dem == 1:
                            pFeature.SetField("elev", dElevation_mean )
                            pFeature.SetField("elev0", dElevation_profile0 )

                        pLayer.CreateFeature(pFeature)

                else:
                    #this hole was added already, but we need to update the neighbor information
                    pCell.aNeighbor_land.append(lCellID)
                    pCell.nNeighbor_land = pCell.nNeighbor_land + 1
                    pCell.aNeighbor_land_virtual = None
                    pCell.nNeighbor_land_virtual = 0

                    pass

        #now update again because some cell has more than one virutal land neighbor, but now none of them is virtual
        #this fix will move virtual land neighbor back to land neighbor
        #the ocean neighbor will remain unchanged
        ncell = len(aMpas_middle)
        for i in range(ncell):
            pCell = aMpas_middle[i]
            aNeighbor_land = pCell.aNeighbor_land           
            aNeighbor_land_virtual_update = list()
            aNeighbor_land_virtual = pCell.aNeighbor_land_virtual
            nNeighbor_land_virtual = pCell.nNeighbor_land_virtual
            nNeighbor_land_update = nNeighbor_land 
            for j in range(nNeighbor_land_virtual):
                lNeighbor = int(aNeighbor_land_virtual[j])
                if lNeighbor in aCellID:
                    #this cell is actually not virtual anymore
                    nNeighbor_land_update = nNeighbor_land_update + 1 
                    aNeighbor_land.append(lNeighbor)
                else:
                    aNeighbor_land_virtual_update.append(lNeighbor)
                    
                

            pCell.aNeighbor_land = aNeighbor_land
            pCell.nNeighbor_land= len(aNeighbor_land)   
            pCell.aNeighbor_land_virtual = aNeighbor_land_virtual_update   
            pCell.nNeighbor_land_virtual = len(aNeighbor_land_virtual_update)
            aMpas_out.append(pCell)


        


    return aMpas_out
