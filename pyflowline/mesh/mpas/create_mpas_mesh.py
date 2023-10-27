import os
import math
import importlib
import numpy as np
from osgeo import ogr, osr, gdal
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
    import netCDF4 as nc

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
        pLayer.CreateField(ogr.FieldDefn('cellid', ogr.OFTInteger64)) #long type for high resolution
        pLayer.CreateField(ogr.FieldDefn('longitude', ogr.OFTReal)) #long type for high resolution
        pLayer.CreateField(ogr.FieldDefn('latitude', ogr.OFTReal)) #long type for high resolution
        pArea_field = ogr.FieldDefn('area', ogr.OFTReal)
        pArea_field.SetWidth(20)
        pArea_field.SetPrecision(2)
        pLayer.CreateField(pArea_field)
        if iFlag_use_mesh_dem == 1:
            pLayer.CreateField(ogr.FieldDefn('elevation_mean', ogr.OFTReal)) #float type for high resolution
            pLayer.CreateField(ogr.FieldDefn('elevation_profile0', ogr.OFTReal)) #float type for high resolution
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
    #aIndexToEdgeID = indexToEdgeID0[:]
    #aIndexToVertexID = indexToVertexID0[:]
    aBed_elevation = bed_elevation0[:]
    aIce_thickness = ice_thickness0[:]
    aCellArea = areaCell0[:]
    aDcEdge = dcEdge0[:]
    aBed_elevation_profile = bed_elevation_profile0[:]  #elevation    
    ncell = len(aIndexToCellID) 
    aMpas = list()
    aMpas_dict = dict()
    lCellIndex=0

    #add a mpas cell into a list
    def add_cell_into_list(aList, i, lCellID, dArea, dElevation_mean, dElevation_profile0, aCoords  ):
        dLon = convert_360_to_180 (aLongitudeCell[i])
        dLat =  (aLatitudeCell[i])        
        if dLon > 100:
            print('Warning: longitude > 100')
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
        if pmpas.nNeighbor != pmpas.nVertex:  #this cell is next to the ocean boundary            
            pmpas.nNeighbor_land = pmpas.nNeighbor 
            pmpas.nNeighbor_ocean = pmpas.nVertex - pmpas.nNeighbor
            pmpas.aNeighbor_land=aNeighborIndex
            pmpas.nNeighbor_land=len(aNeighborIndex)
            print(lCellID)
            print('Warning: nNeighbor != nVertex')
        else: #this cell is not at the the land-ocean mask coastal line
            pmpas.nNeighbor_land = pmpas.nNeighbor 
            pmpas.nNeighbor_ocean = 0
            pmpas.aNeighbor_land = aNeighborIndex        
            pmpas.nNeighbor_land=len(aNeighborIndex)

        aDistance=list()
        for j in range(pmpas.nNeighbor):
            #find shared edge
            lEdgeID= aEdgeIndex[j]
            lIndex = lEdgeID-1
            dDistance = aDcEdge[lIndex]
            aDistance.append(dDistance)
            pass

        #this contains all the original mpas neighbor distance
        pmpas.aNeighbor_distance = aDistance 
        aList.append(pmpas)
        return aList

    if iFlag_antarctic == 1:
        iFlag_remove_ice = 0
        #if it is antarctic, we dont need the boundary
        for i in range(ncell):
            #center
            dLon = convert_360_to_180 (aLongitudeCell[i])
            dLat =  (aLatitudeCell[i])
            
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
                aMpas_dict[lCellID] = lCellIndex
                lCellIndex = lCellIndex + 1
                #save mesh cell
                if iFlag_save_mesh_in ==1:                
                    pFeature.SetGeometry(pPolygon)
                    pFeature.SetField("cellid", int(lCellID) )
                    pFeature.SetField("longitude", dLon )
                    pFeature.SetField("latitude", dLat )
                    pFeature.SetField("area", dArea )
                    if iFlag_use_mesh_dem == 1:
                        pFeature.SetField("elevation_mean", dElevation_mean )
                        pFeature.SetField("elevation_profile0", dElevation_profile0 )

                    pLayer.CreateFeature(pFeature)
            
    else:
        iFlag_remove_ice = 1
        for i in range(ncell):
            #center
                  
            #vertex
            aVertexOnCellIndex = np.array(aVertexOnCell[i,:])
            dummy0 = np.where(aVertexOnCellIndex > 0)
            aVertexIndex = aVertexOnCellIndex[dummy0]
            aLonVertex = aLongitudeVertex[aVertexIndex-1]
            aLatVertex = aLatitudeVertex[aVertexIndex-1]
            nVertex = len(aLonVertex)
            #first check if it is within the boundary
            
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
                dLon_min = np.min(aCoords[:,0])
                dLon_max = np.max(aCoords[:,0])
                if np.abs(dLon_min-dLon_max) > 100: #this polygon cross international date line
                    #print('Warning: longitude > 100')
                    pass
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
                aMpas_dict[lCellID] = lCellIndex
                lCellIndex = lCellIndex + 1
                #save mesh cell
                if iFlag_save_mesh_in == 1:      
                    dLon = convert_360_to_180 (aLongitudeCell[i])
                    dLat =  (aLatitudeCell[i])            
                    pFeature.SetGeometry(pPolygon)
                    pFeature.SetField("cellid", int(lCellID) )                   
                    pFeature.SetField("longitude", dLon )
                    pFeature.SetField("latitude", dLat )
                    pFeature.SetField("area", dArea )
                    if iFlag_use_mesh_dem == 1:
                        pFeature.SetField("elevation_mean", dElevation_mean )
                        pFeature.SetField("elevation_profile0", dElevation_profile0 )

                    pLayer.CreateFeature(pFeature)


    #for maps we need to clean some cell because they were not actually in the domain
    #besides, we need to add some smal holes back
    #to do this, we need two steps.
    

    if iFlag_global_in == 1:
        aMpas_out = aMpas
    else:
        iFlag_fill_hole = 0        
        aMpas_out = list()
        ncell = len(aMpas)
        #generate the list of cell ID that are already certain        
        if iFlag_fill_hole == 1:  
            
            #first update neighbor information because some cell should have vitual land neighbor (not present in the mesh)
            #this operation does not increase the number of cells, but it update the neighbor information
            #specifically, it divided the land neighbor into two parts: land and virtual land         
            for pCell in aMpas:               
                aNeighbor_land = pCell.aNeighbor_land   #including both holes and maps land cutoff by boundary
                aNeighbor_land_update = list()
                aNeighbor_land_virtual = list()
                for lNeighbor in aNeighbor_land: #loop all land neighbors                    
                    if lNeighbor in aMpas_dict:                        
                        aNeighbor_land_update.append(lNeighbor)
                    else:
                        #a hole or boundary mpas land cell
                        aNeighbor_land_virtual.append(lNeighbor)

                pCell.aNeighbor_land = aNeighbor_land_update
                pCell.nNeighbor_land= len(aNeighbor_land_update)   
                pCell.aNeighbor_land_virtual = aNeighbor_land_virtual   
                pCell.nNeighbor_land_virtual = len(aNeighbor_land_virtual)
                pass
            

                #distance remains unchanged since we just have missing cells.

            #now add back small holes 
            #this operation will increase the number of cells  
            #it will also update the neighbor information for some cells, 
            #not all cells will be updated (because some cells have 2+ virtual land neighbors)
            for pCell in aMpas:                     
                if pCell.nNeighbor_land_virtual == 1:  #only one virtual land means it is likely next to a hole 
                    lNeighbor_hole = pCell.aNeighbor_land_virtual[0]
                    j = lNeighbor_hole-1
                    dLon = convert_360_to_180 (aLongitudeCell[j])
                    dLat =  aLatitudeCell[j]                    
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

                    x1 = convert_360_to_180(aLonVertex[0])
                    y1 = aLatVertex[0]
                    ring.AddPoint(x1, y1) #double check            
                    pPolygon = ogr.Geometry(ogr.wkbPolygon)
                    pPolygon.AddGeometry(ring)

                    lCellID = int(aIndexToCellID[j])
                    dElevation_mean = float(aBed_elevation[j])
                    dElevation_profile0 = float(aBed_elevation_profile[j,0])
                    dArea = float(aCellArea[j])

                    if lCellID not in aMpas_dict:   
                        aMpas = add_cell_into_list(aMpas, j, lCellID, dArea, dElevation_mean, dElevation_profile0, aCoords )
                        aMpas_dict[lCellID] = lCellIndex
                        lCellIndex = lCellIndex + 1
                        #now we need to update the neightboring information as well
                        pCell.aNeighbor_land.append(lCellID)
                        pCell.nNeighbor_land = pCell.nNeighbor_land + 1
                        pCell.aNeighbor_land_virtual = None
                        pCell.nNeighbor_land_virtual = 0
                        #save mesh cell
                        if iFlag_save_mesh_in ==1:                
                            pFeature.SetGeometry(pPolygon)
                            pFeature.SetField("cellid", int(lCellID) )
                            pFeature.SetField("longitude", dLon )
                            pFeature.SetField("latitude", dLat )
                            pFeature.SetField("area", dArea )
                            if iFlag_use_mesh_dem == 1:
                                pFeature.SetField("elevation_mean", dElevation_mean )
                                pFeature.SetField("elevation_profile0", dElevation_profile0 )

                            pLayer.CreateFeature(pFeature)

                    else:
                        #this hole was added already, but we need to update the neighbor information
                        pCell.aNeighbor_land.append(lCellID)
                        pCell.nNeighbor_land = pCell.nNeighbor_land + 1
                        pCell.aNeighbor_land_virtual = None
                        pCell.nNeighbor_land_virtual = 0
                        pass

                #how about distance? still unchanged, but the orders are changed

            #now update again because some cell has more than one virutal land neighbor, but now none of them is virtual anymore
            #this fix will move virtual land neighbor back to land neighbor
            #the ocean neighbor will remain unchanged
        
            for pCell in aMpas:            
                aNeighbor_land_update = list()   
                aNeighbor_land = pCell.aNeighbor_land  
                aNeighbor_land_virtual_update = list()      
                aNeighbor_land_virtual = pCell.aNeighbor_land_virtual                        
                for lNeighbor in aNeighbor_land:     
                    if lNeighbor in aMpas_dict:
                        aNeighbor_land_update.append(lNeighbor)                    
                        pass
                    else:
                        #this is a land cell in mpas, but it may be clipped by boundary
                        pass

                #for book keeping only        
                for lNeighbor in aNeighbor_land_virtual:         
                    if lNeighbor in aMpas_dict:
                        #this cell is actually not virtual anymore                    
                        aNeighbor_land_update.append(lNeighbor)
                    else:
                        aNeighbor_land_virtual_update.append(lNeighbor)
                        pass

                pCell.aNeighbor_land = aNeighbor_land_update
                pCell.nNeighbor_land= len(aNeighbor_land_update)   
                pCell.aNeighbor_land_virtual = aNeighbor_land_virtual_update   #for book keeping only
                pCell.nNeighbor_land_virtual = len(aNeighbor_land_virtual_update)
                aMpas_out.append(pCell)
        else:
            #no hole filling applied
            #still need to get rid cell that are not in the domain
            for pCell in aMpas:       
                aNeighbor = pCell.aNeighbor         
                aNeighbor_land_update = list() 
                for lNeighbor in aNeighbor:              
                    if lNeighbor in aMpas_dict:           
                        aNeighbor_land_update.append(lNeighbor)

                #for latlon, there is no ocean concept
                pCell.aNeighbor = aNeighbor_land_update
                pCell.nNeighbor= len(aNeighbor_land_update)     
                pCell.aNeighbor_land = aNeighbor_land_update       
                pCell.nNeighbor_land= len(aNeighbor_land_update)            
                pCell.nNeighbor_ocean = pCell.nVertex - pCell.nNeighbor_land
                aMpas_out.append(pCell)
          
            pass

    return aMpas_out
