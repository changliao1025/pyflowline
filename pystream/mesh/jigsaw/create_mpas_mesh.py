import os, sys
import math
import numpy as np
from netCDF4 import Dataset

from osgeo import ogr, osr, gdal, gdalconst

from pystream.shared.mpas import pympas
from shapely.wkt import loads
from pyearth.gis.location.convert_lat_lon_range import convert_180_to_360,convert_360_to_180

from pystream.format.convert_coordinates_to_cell import convert_coordinates_to_cell
from pystream.format.convert_attribute_to_cell import convert_attribute_to_cell

def create_mpas_mesh(sFilename_mesh_netcdf, dLatitude_top, dLatitude_bot, dLongitude_left, dLongitude_right,sFilename_mesh):
    
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
    #pDataset_shapefile = pDriver_shapefile.Open(sFilename_shapefile, 0)
    #pLayer_shapefile = pDataset_shapefile.GetLayer(0)
    #pSrs = pLayer_shapefile.GetSpatialRef()
    pSrs = osr.SpatialReference()  
    pSrs.ImportFromEPSG(4326)    # WGS84 lat/lon

    
    #geojson
    pDataset = pDriver_shapefile.CreateDataSource(sFilename_mesh)

    pLayer = pDataset.CreateLayer('cell', pSrs, ogr.wkbPolygon)
    # Add one attribute
    pLayer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger64)) #long type for high resolution
    
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
    ncell = len(aIndexToCellID)
    
 
    aMpas = list()
    for i in range(ncell):
        dLat = convert_360_to_180 (aLatitudeCell[i])
        dLon = convert_360_to_180 (aLongitudeCell[i])

        if dLat > dLatitude_bot and dLat < dLatitude_top and dLon > dLongitude_left and dLon < dLongitude_right:


            #get cell edge
            lCellID = int(aIndexToCellID[i])
            aCellOnCellIndex = np.array(aCellsOnCell[i,:])
            aEdgesOnCellIndex = np.array(aEdgesOnCell[i,:])
            aVertexOnCellIndex = np.array(aVertexOnCell[i,:])

            dummy0 = np.where(aVertexOnCellIndex > 0)
            aVertexIndex = aVertexOnCellIndex[dummy0]
            aEdgeIndex= aEdgesOnCellIndex[dummy0]
            aNeighborIndex= aCellOnCellIndex[dummy0]

            #
            aVertexIndexOnEdge = np.array(aVertexOnEdge0[aEdgeIndex-1,:]).astype((np.int))

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
            x1 = convert_360_to_180(aLonVertex[0])
            y1 = aLatVertex[0]
            ring.AddPoint(x1, y1) #double check            

            pPolygon = ogr.Geometry(ogr.wkbPolygon)
            pPolygon.AddGeometry(ring)

            pFeature.SetGeometry(pPolygon)
            pFeature.SetField("id", int(lCellID) )
            pLayer.CreateFeature(pFeature)
            
            pmpas = convert_attribute_to_cell(4, aVertexIndex, aEdgeIndex, aVertexIndexOnEdge ,aCoords)
           
            pmpas.lCellID = lCellID
            
            pmpas.aNeighbor=aNeighborIndex
            pmpas.nNeighbor=nVertex
            aMpas.append(pmpas)
    
            #get vertex

        pass

    return aMpas
