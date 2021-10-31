import os, sys
import math
import numpy as np
from netCDF4 import Dataset

from osgeo import ogr, osr, gdal, gdalconst

from pyflowline.shared.mpas import pympas
from shapely.wkt import loads
from pyearth.gis.location.convert_lat_lon_range import convert_180_to_360,convert_360_to_180

from pyflowline.format.convert_coordinates_to_cell import convert_pcs_coordinates_to_cell
from pyflowline.format.convert_attribute_to_cell import convert_gcs_attribute_to_cell

def create_mpas_mesh(iFlag_use_mesh_dem, iFlag_save_mesh, \
    dLatitude_top, dLatitude_bot, dLongitude_left, dLongitude_right,\
     sFilename_mesh_netcdf, sFilename_mesh):
    
    if (os.path.exists(sFilename_mesh_netcdf)):
        pass
    else:
        print('Mesh file does not exist!')
        exit
    
    if os.path.exists(sFilename_mesh): 
        #delete it if it exists
        os.remove(sFilename_mesh)

    pDatasets_in = Dataset(sFilename_mesh_netcdf)

    netcdf_format = pDatasets_in.file_format
    pDriver_geojson = ogr.GetDriverByName('GeoJSON')
    pDriver_shapefile = ogr.GetDriverByName('ESRI Shapefile')
  
    pSpatialRef_gcs = osr.SpatialReference()  
    pSpatialRef_gcs.ImportFromEPSG(4326)    # WGS84 lat/lon

    
    #geojson
    if iFlag_save_mesh ==1:
        pDataset = pDriver_shapefile.CreateDataSource(sFilename_mesh)

        pLayer = pDataset.CreateLayer('cell', pSpatialRef_gcs, ogr.wkbPolygon)
        # Add one attribute
        pLayer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger64)) #long type for high resolution
        pLayer.CreateField(ogr.FieldDefn('lon', ogr.OFTReal)) #long type for high resolution
        pLayer.CreateField(ogr.FieldDefn('lat', ogr.OFTReal)) #long type for high resolution

        if iFlag_use_mesh_dem == 1:
            pLayer.CreateField(ogr.FieldDefn('elev', ogr.OFTReal)) #float type for high resolution
        else:

            pass
    
        pLayerDefn = pLayer.GetLayerDefn()
        pFeature = ogr.Feature(pLayerDefn)
        

    #write new netcdf
    for sKey, aValue in pDatasets_in.variables.items():        
        print(aValue.datatype)
        print(aValue.dimensions)

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

        if sKey == 'ice_elevation':
            ice_elevation0 = aValue 
        else:
            pass

        if sKey == 'areaCell':
            areaCell0 = aValue 
        else:
            pass
        

    aLatitudeVertex = latVertex0[:] / math.pi * 180
    aLongitudeVertex = lonVertex0[:] / math.pi * 180
    
    #conver unit 
    aLatitudeCell = latCell0[:] / math.pi * 180
    aLongitudeCell = lonCell0[:] / math.pi * 180

    aCellsOnCell = cellsOnCell0[:]

    aEdgesOnCell= edgesOnCell0[:]
    aVertexOnCell = verticesOnCell0[:]
    aVertexOnEdge0 = verticesOnEdge0[:]
    
    aIndexToCellID = indexToCellID0[:]
    aIndexToEdgeID = indexToEdgeID0[:]
    aIndexToVertexID = indexToVertexID0[:]

    aBed_elevation = bed_elevation0[:]
    aIce_thickness = ice_elevation0[:]
    aCellArea = areaCell0[:]
    
    ncell = len(aIndexToCellID)
 
    aMpas = list()
    for i in range(ncell):
        dLat = convert_360_to_180 (aLatitudeCell[i])
        dLon = convert_360_to_180 (aLongitudeCell[i])

        if dLat > dLatitude_bot and dLat < dLatitude_top and dLon > dLongitude_left and dLon < dLongitude_right:

            #get cell edge
            lCellID = int(aIndexToCellID[i])
            dElevation = float(aBed_elevation[i])
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

                if iFlag_save_mesh ==1:
                    x1 = convert_360_to_180(aLonVertex[0])
                    y1 = aLatVertex[0]
                    ring.AddPoint(x1, y1) #double check            
                    pPolygon = ogr.Geometry(ogr.wkbPolygon)
                    pPolygon.AddGeometry(ring)
                    pFeature.SetGeometry(pPolygon)
                    pFeature.SetField("id", int(lCellID) )
                    pFeature.SetField("lon", dLon )
                    pFeature.SetField("lat", dLat )
                    if iFlag_use_mesh_dem == 1:
                        pFeature.SetField("elev", dElevation )

                    pLayer.CreateFeature(pFeature)

                pmpas = convert_gcs_attribute_to_cell(4, dLon, dLat, aCoords, aVertexIndex, aEdgeIndex, aVertexIndexOnEdge)
                #calculate area 
                #pmpas.calculate_cell_area()
                pmpas.dArea = dArea
                pmpas.calculate_edge_length()
                pmpas.lCellID = lCellID
                pmpas.dElevation  = dElevation
                pmpas.aNeighbor=aNeighborIndex
                pmpas.nNeighbor=len(aNeighborIndex)
                pmpas.aNeighbor_land=aNeighborIndex
                pmpas.nNeighbor_land=len(aNeighborIndex)
                aMpas.append(pmpas)
                #get vertex

        pass

    #for maps we need to clean some cell because they were not actually in the domain
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
