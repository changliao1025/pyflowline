
import numpy as np
from osgeo import ogr, gdal
from shapely.wkt import loads

from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_cell

def read_mesh_json(iMesh_type_in, sFilename_mesh_in):
    """
    convert a shpefile to json format.
    This function should be used for stream flowline only.
    """

    aCell_out=list()
    pDriver_json = ogr.GetDriverByName('GeoJSON') 
   
    pDataset_mesh = pDriver_json.Open(sFilename_mesh_in, gdal.GA_ReadOnly)
    pLayer_mesh = pDataset_mesh.GetLayer(0)
    pSpatial_reference_out = pLayer_mesh.GetSpatialRef()

    ldefn = pLayer_mesh.GetLayerDefn()
    schema =list()
    for n in range(ldefn.GetFieldCount()):
        fdefn = ldefn.GetFieldDefn(n)
        schema.append(fdefn.name)
    if 'iseg' in schema:
        iFlag_segment = 1
    else:
        iFlag_segment = 0   

    #we also need to spatial reference
    for pFeature_mesh in pLayer_mesh:
       
        #pFeature_mesh= pLayer_mesh.GetFeature(i)
        pGeometry_mesh = pFeature_mesh.GetGeometryRef()        
        dummy0 = loads( pGeometry_mesh.ExportToWkt() )
        aCoords_gcs = dummy0.exterior.coords
        aCoords_gcs= np.array(aCoords_gcs)       

        lCellID = pFeature_mesh.GetField("id")
        dLon = pFeature_mesh.GetField("lon")
        dLat = pFeature_mesh.GetField("lat")
        
        dArea = pFeature_mesh.GetField("area")
        if iMesh_type_in == 4:
            dElevation_mean = pFeature_mesh.GetField("elev")
            dElevation_profile0 = pFeature_mesh.GetField("elev0")
        #convert geometry to edge
        pGeometrytype_mesh = pGeometry_mesh.GetGeometryName()
        if(pGeometrytype_mesh == 'POLYGON'):            
            pCell = convert_gcs_coordinates_to_cell(iMesh_type_in, dLon, dLat, aCoords_gcs)     
            pCell.lCellID = lCellID #this information is saved in shapefile            
            pCell.dArea = dArea #pCell.calculate_cell_area()
            pCell.dLength = pCell.calculate_edge_length()
            pCell.dLength_flowline = pCell.dLength
            if iMesh_type_in == 4:
                pCell.dElevation_mean = dElevation_mean
                pCell.dElevation_profile0 = dElevation_profile0

            aCell_out.append(pCell)

    return aCell_out, pSpatial_reference_out



def read_mesh_shapefile(sFilename_mesh_in):
    """
    convert a shpefile to json format.
    This function should be used for stream flowline only.
    """
    aMesh_out=list()

    
    pDriver_shapefile = ogr.GetDriverByName('ESRI Shapefile')
   
    pDataset_shapefile = pDriver_shapefile.Open(sFilename_mesh_in, gdal.GA_ReadOnly)
    pLayer_shapefile = pDataset_shapefile.GetLayer(0)
    pSpatial_reference_out = pLayer_shapefile.GetSpatialRef()

    #we also need to spatial reference

    return aMesh_out, pSpatial_reference_out