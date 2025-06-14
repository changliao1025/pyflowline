import os
from osgeo import ogr, osr
from pyearth.gis.geometry.calculate_polygon_area  import calculate_polygon_area

def export_vertex_to_geojson(aVertex_in,
        sFilename_json_in,
        iFlag_projected_in=None,
        pSpatial_reference_in=None,
        aAttribute_data=None):
    """
    Convert a shapefile to json format.
    This function should be used for stream flowline only.

    Args:
        aVertex_in (_type_): _description_
        sFilename_json_in (_type_): _description_
        iFlag_projected_in (_type_, optional): _description_. Defaults to None.
        pSpatial_reference_in (_type_, optional): _description_. Defaults to None.
        aAttribute_data (_type_, optional): _description_. Defaults to None.
    """


    if os.path.exists(sFilename_json_in):
        os.remove(sFilename_json_in)

    iFlag_projected_in = 0 if iFlag_projected_in is None else 1

    if  pSpatial_reference_in is None:
        pSpatial_reference_in = osr.SpatialReference()
        pSpatial_reference_in.ImportFromEPSG(4326)    # WGS84 lat/lon


    iFlag_attribute = 1 if aAttribute_data is not None else 0
    aAttribute = aAttribute_data if aAttribute_data is not None else []

    #nVertex = len(aVertex_in)
    pDriver = ogr.GetDriverByName('GeoJSON')
    pDataset_json = pDriver.CreateDataSource(sFilename_json_in)
    pLayer_json = pDataset_json.CreateLayer('vertex', pSpatial_reference_in, ogr.wkbPoint)
    # Add one attribute
    pLayer_json.CreateField(ogr.FieldDefn('pointid', ogr.OFTInteger64)) #long type for high resolution
    if iFlag_attribute ==1:
        pLayer_json.CreateField(ogr.FieldDefn('connectivity', ogr.OFTInteger64)) #long type for high resolution
        pass

    pLayerDefn = pLayer_json.GetLayerDefn()
    pFeature_out = ogr.Feature(pLayerDefn)
    #lID = 0
    #for i in range(nVertex):
    #    pVertex = aVertex_in[i]
    #    pPoint = ogr.Geometry(ogr.wkbPoint)
    #    if iFlag_projected_in ==1:
    #        pPoint.AddPoint(pVertex.dx, pVertex.dy)
    #        pass
    #    else:
    #        pPoint.AddPoint(pVertex.dLongitude_degree, pVertex.dLatitude_degree)
    #        pass
    #    pGeometry_out = ogr.CreateGeometryFromWkb(pPoint.ExportToWkb())
    #    pFeature_out.SetGeometry(pGeometry_out)
    #    pFeature_out.SetField("pointid", lID)
    #    if iFlag_attribute ==1:
    #        pFeature_out.SetField("connectivity", int(aAttribute[i]) )
    #
    #    pLayer_json.CreateFeature(pFeature_out)
    #    lID =  lID + 1
    #    pass

    for lID, pVertex in enumerate(aVertex_in):
        pPoint = ogr.Geometry(ogr.wkbPoint)
        if iFlag_projected_in == 1:
            pPoint.AddPoint(pVertex.dx, pVertex.dy)
        else:
            pPoint.AddPoint(pVertex.dLongitude_degree, pVertex.dLatitude_degree)

        pGeometry_out = ogr.CreateGeometryFromWkb(pPoint.ExportToWkb())
        pFeature_out.SetGeometry(pGeometry_out)
        pFeature_out.SetField("pointid", lID + 1)

        if iFlag_attribute == 1:
            pFeature_out.SetField("connectivity", int(aAttribute[lID]))

        pLayer_json.CreateFeature(pFeature_out)

    pDataset_json.FlushCache()
    pDataset_json = pLayer_json = pFeature_out  = None

    return

def export_vertex_to_shapefile(aVertex_in, sFilename_shapefile_out,\
    pSpatial_reference_in,
    iFlag_projected_in,
    aAttribute_data=None):
    """
    convert a shpefile to json format.
    This function should be used for stream flowline only.
    """

    if os.path.exists(sFilename_shapefile_out):
        os.remove(sFilename_shapefile_out)
        pass

    if aAttribute_data is not None:
        aAttribute= aAttribute_data
        iFlag_attribute =1
    else:
        iFlag_attribute=0

    nVertex = len(aVertex_in)
    pDriver = ogr.GetDriverByName('ESRI Shapefile')
    pDataset_json = pDriver.CreateDataSource(sFilename_shapefile_out)
    pLayer_json = pDataset_json.CreateLayer('vertex', pSpatial_reference_in, ogr.wkbPoint)
    # Add one attribute
    pLayer_json.CreateField(ogr.FieldDefn('id', ogr.OFTInteger64)) #long type for high resolution
    if iFlag_attribute ==1:
        pLayer_json.CreateField(ogr.FieldDefn('con', ogr.OFTInteger64)) #long type for high resolution
        pass

    pLayerDefn = pLayer_json.GetLayerDefn()
    pFeature_out = ogr.Feature(pLayerDefn)
    lID = 0
    for i in range(nVertex):
        pVertex = aVertex_in[i]
        pPoint = ogr.Geometry(ogr.wkbPoint)
        if iFlag_projected_in ==1:
            #dummy1= Point( pVertex.dx, pVertex.dy )
            pPoint.AddPoint(pVertex.dx, pVertex.dy)
            pass
        else:
            #dummy1= Point( pVertex.dLongitude_degree, pVertex.dLatitude_degree )
            pPoint.AddPoint(pVertex.dLongitude_degree, pVertex.dLatitude_degree)
            pass

        #pGeometry_out = ogr.CreateGeometryFromWkb(dummy1.wkb)
        pGeometry_out = ogr.CreateGeometryFromWkb(pPoint.ExportToWkb())
        pFeature_out.SetGeometry(pGeometry_out)
        pFeature_out.SetField("id", lID)
        if iFlag_attribute ==1:
            pFeature_out.SetField("con", int(aAttribute[i]) )

        # Add new pFeature_shapefile to output Layer
        pLayer_json.CreateFeature(pFeature_out)
        lID =  lID + 1
        pass

    pDataset_json.FlushCache()
    pDataset_json = pLayer_json = pFeature_out  = None

    return

def export_vertex_as_polygon(aVertex_in, sFilename_out):
    pDriver = ogr.GetDriverByName('GEOJSON')
    if os.path.exists(sFilename_out):
        os.remove(sFilename_out)

    pDataset = pDriver.CreateDataSource(sFilename_out)
    pSrcSpatialRef = osr.SpatialReference()
    pSrcSpatialRef.ImportFromEPSG(4326)
    pLayer = pDataset.CreateLayer('buffer_ploygon', pSrcSpatialRef, geom_type=ogr.wkbPolygon)
    aLon = list()
    aLat = list()
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for pVertex in aVertex_in:
        ring.AddPoint(pVertex.dLongitude_degree, pVertex.dLatitude_degree)
        aLon.append(pVertex.dLongitude_degree)
        aLat.append(pVertex.dLatitude_degree)

    ring.AddPoint(aVertex_in[0].dLongitude_degree, aVertex_in[0].dLatitude_degree)
    pArea_field = ogr.FieldDefn('area', ogr.OFTReal)
    pArea_field.SetWidth(20)
    pArea_field.SetPrecision(2)
    pLayer.CreateField(pArea_field)
    pLayerDefn = pLayer.GetLayerDefn()
    pFeature = ogr.Feature(pLayerDefn)
    #add the first point to close the ring
    dArea = calculate_polygon_area(aLon, aLat)
    pPolygon = ogr.Geometry(ogr.wkbPolygon)
    pPolygon.AddGeometry(ring)
    pFeature.SetGeometry(pPolygon)
    pFeature.SetField("area", dArea )
    pLayer.CreateFeature(pFeature)
    #flush the cache and close the file
    pDataset.FlushCache()
    pDataset.Destroy()
    pSrcSpatialRef = None


    return



