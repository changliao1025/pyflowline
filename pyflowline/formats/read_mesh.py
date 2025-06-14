import os
import numpy as np
from osgeo import ogr, gdal
#from shapely.wkt import loads
from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_cell
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates
def read_mesh_json(iMesh_type_in, sFilename_mesh_in):
    """
    convert a shpefile to json format.
    This function should be used for stream flowline only.
    """
    iReturn_code = 1
    if os.path.isfile(sFilename_mesh_in):
        pass
    else:
        print('This mesh file does not exist: ', sFilename_mesh_in )
        iReturn_code = 0
        return iReturn_code

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
    if 'stream_segment' in schema:
        iFlag_segment = 1
    else:
        iFlag_segment = 0

    #we also need to spatial reference
    for pFeature_mesh in pLayer_mesh:
        pGeometry_mesh = pFeature_mesh.GetGeometryRef()
        #dummy0 = loads( pGeometry_mesh.ExportToWkt() )
        #aCoords_gcs = dummy0.exterior.coords
        #aCoords_gcs= np.array(aCoords_gcs)
        aCoords_gcs = get_geometry_coordinates(pGeometry_mesh)
        lCellID = pFeature_mesh.GetField("cellid")
        dLongitude_center = pFeature_mesh.GetField("longitude")
        dLatitude_center = pFeature_mesh.GetField("latitude")
        dArea = pFeature_mesh.GetField("area")
        if iMesh_type_in == 4:
            dElevation_mean = pFeature_mesh.GetField("elevation_mean")
            dElevation_profile0 = pFeature_mesh.GetField("elevation_profile0")

        pGeometrytype_mesh = pGeometry_mesh.GetGeometryName()
        if(pGeometrytype_mesh == 'POLYGON'):
            pCell = convert_gcs_coordinates_to_cell(iMesh_type_in, dLongitude_center, dLatitude_center, aCoords_gcs)
            pCell.lCellID = lCellID
            pCell.dArea = dArea
            pCell.dLength = pCell.calculate_edge_length()
            pCell.dLength_flowline = pCell.dLength
            if iMesh_type_in == 4:
                pCell.dElevation_mean = dElevation_mean
                pCell.dElevation_profile0 = dElevation_profile0

            aCell_out.append(pCell)

    return aCell_out, pSpatial_reference_out

def read_mesh_json_w_topology(iMesh_type_in, sFilename_mesh_in):
    """

    This function should be used for stream flowline only.
    """
    iReturn_code = 1
    if os.path.isfile(sFilename_mesh_in):
        pass
    else:
        print('This mesh file does not exist: ', sFilename_mesh_in )
        iReturn_code = 0
        return iReturn_code

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
    if 'stream_segment' in schema:
        iFlag_segment = 1
    else:
        iFlag_segment = 0

    #we also need to spatial reference
    for pFeature_mesh in pLayer_mesh:
        pGeometry_mesh = pFeature_mesh.GetGeometryRef()
        #dummy0 = loads( pGeometry_mesh.ExportToWkt() )
        #aCoords_gcs = dummy0.exterior.coords
        #aCoords_gcs= np.array(aCoords_gcs)
        aCoords_gcs = get_geometry_coordinates(pGeometry_mesh)
        lCellID = pFeature_mesh.GetField("cellid")
        dLongitude_center = pFeature_mesh.GetField("longitude")
        dLatitude_center = pFeature_mesh.GetField("latitude")
        dArea = pFeature_mesh.GetField("area")
        if iMesh_type_in == 4:
            dElevation_mean = pFeature_mesh.GetField("elevation_mean")
            dElevation_profile0 = pFeature_mesh.GetField("elevation_profile0")

        pGeometrytype_mesh = pGeometry_mesh.GetGeometryName()
        if(pGeometrytype_mesh == 'POLYGON'):
            pCell = convert_gcs_coordinates_to_cell(iMesh_type_in, dLongitude_center, dLatitude_center, aCoords_gcs)
            pCell.lCellID = lCellID
            pCell.dArea = dArea
            pCell.dLength = pCell.calculate_edge_length()
            pCell.dLength_flowline = pCell.dLength
            if iMesh_type_in == 4:
                pCell.dElevation_mean = dElevation_mean
                pCell.dElevation_profile0 = dElevation_profile0

            aCell_out.append(pCell)

    return aCell_out, pSpatial_reference_out