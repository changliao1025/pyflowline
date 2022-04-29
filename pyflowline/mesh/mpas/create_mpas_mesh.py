import os
import math
import numpy as np
from osgeo import ogr, osr
from netCDF4 import Dataset
from pyflowline.classes.mpas import pympas
from pyflowline.formats.convert_attributes import convert_gcs_attributes_to_cell
from pyflowline.algorithms.auxiliary.gdal_functions import convert_360_to_180

def create_mpas_mesh(iFlag_global_in, \
    iFlag_use_mesh_dem, \
        iFlag_save_mesh_in, \
        dLongitude_left_in, dLongitude_right_in,\
    dLatitude_top_in, dLatitude_bot_in, \
     sFilename_mesh_netcdf_in, \
         sFilename_output_in):
    
    if (os.path.exists(sFilename_mesh_netcdf_in)):
        pass
    else:
        print('Mesh file does not exist!')
        exit
    
    if os.path.exists(sFilename_output_in):  
        os.remove(sFilename_output_in)

    pDatasets_in = Dataset(sFilename_mesh_netcdf_in)

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

    #write new netcdf
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
    for i in range(ncell):
        dLat = convert_360_to_180 (aLatitudeCell[i])
        dLon = convert_360_to_180 (aLongitudeCell[i])
        if dLat > dLatitude_bot_in and dLat < dLatitude_top_in and dLon > dLongitude_left_in and dLon < dLongitude_right_in:
            #get cell edge
            lCellID = int(aIndexToCellID[i])
            dElevation_mean = float(aBed_elevation[i])
            dElevation_profile0 = float(aBed_elevation_profile[i,0])
            dThickness_ice = float( aIce_thickness[i] )
            dArea = float(aCellArea[i])
            if dThickness_ice > 0:
                continue
            else:
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
                aLonVertex = aLongitudeVertex[aVertexIndex-1]
                aLatVertex = aLatitudeVertex[aVertexIndex-1]
                nVertex = len(aLonVertex)
                ring = ogr.Geometry(ogr.wkbLinearRing)
                aCoords = np.full((nVertex,2), -9999.0, dtype=float)

                for j in range(nVertex):
                    x1 = convert_360_to_180(aLonVertex[j])
                    y1 = aLatVertex[j]      
                    ring.AddPoint(x1, y1)
                    aCoords[j,0] = x1
                    aCoords[j,1] = y1
                    pass

                if iFlag_save_mesh_in ==1:
                    x1 = convert_360_to_180(aLonVertex[0])
                    y1 = aLatVertex[0]
                    ring.AddPoint(x1, y1) #double check            
                    pPolygon = ogr.Geometry(ogr.wkbPolygon)
                    pPolygon.AddGeometry(ring)
                    pFeature.SetGeometry(pPolygon)
                    pFeature.SetField("id", int(lCellID) )
                    pFeature.SetField("lon", dLon )
                    pFeature.SetField("lat", dLat )
                    pFeature.SetField("area", dArea )
                    if iFlag_use_mesh_dem == 1:
                        pFeature.SetField("elev", dElevation_mean )
                        pFeature.SetField("elev0", dElevation_profile0 )

                    pLayer.CreateFeature(pFeature)

                pmpas = convert_gcs_attributes_to_cell(4, dLon, dLat, aCoords, aVertexIndex, aEdgeIndex, aVertexIndexOnEdge)               
                pmpas.dArea = dArea
                pmpas.calculate_edge_length()
                pmpas.dLength_flowline = pmpas.dLength_edge #Default
                pmpas.lCellID = lCellID
                pmpas.dElevation_mean  = dElevation_mean
                pmpas.dElevation_profile0 = dElevation_profile0
                pmpas.aNeighbor=aNeighborIndex
                pmpas.nNeighbor=len(aNeighborIndex)
                pmpas.aNeighbor_land=aNeighborIndex
                pmpas.nNeighbor_land=len(aNeighborIndex)
                aDistance=list()
                for i in range(pmpas.nNeighbor):
                    lNeighborID = pmpas.aNeighbor[i]
                    #find shared edge
                    lEdgeID= aEdgeIndex[i]                    
                    lIndex = aIndexToEdgeID[lEdgeID-1]
                    dDistance = aDcEdge[lIndex]
                    aDistance.append(dDistance)
                    pass

                pmpas.aNeighbor_distance = aDistance
                aMpas.append(pmpas)
                #get vertex

        pass

    #for maps we need to clean some cell because they were not actually in the domain

    if iFlag_global_in == 1:
        aMpas_out = aMpas
    else:
        aMpas_out = list()
        ncell = len(aMpas)
        aCellID  = list()
        for i in range(ncell):
            pCell = aMpas[i]
            lCellID = pCell.lCellID
            aCellID.append(lCellID)

        for i in range(ncell):
            pCell = aMpas[i]
            aNeighbor = pCell.aNeighbor
            nNeighbor = pCell.nNeighbor
            aNeighbor_new = list()
            nNeighbor_new = 0 
            for j in range(nNeighbor):
                lNeighbor = int(aNeighbor[j])
                if lNeighbor in aCellID:
                    nNeighbor_new = nNeighbor_new + 1 
                    aNeighbor_new.append(lNeighbor)

            pCell.nNeighbor_land= len(aNeighbor_new)
            pCell.aNeighbor_land = aNeighbor_new
            pCell.nNeighbor_ocean = pCell.nVertex - pCell.nNeighbor_land
            aMpas_out.append(pCell)


    return aMpas_out
