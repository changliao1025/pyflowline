import os
from osgeo import ogr, osr
from shapely.geometry import Point

def export_vertex_to_json(aVertex_in, \
        sFilename_json_in,\
        iFlag_projected_in=None,\
        pSpatial_reference_in=None, \
        aAttribute_data=None):
    """
    convert a shapefile to json format.
    This function should be used for stream flowline only.
    """

    if os.path.exists(sFilename_json_in): 
        os.remove(sFilename_json_in)
        pass
    
    if iFlag_projected_in is None:
        iFlag_projected_in = 0
    else:
        iFlag_projected_in = 1

    if  pSpatial_reference_in is None:        
        pSpatial_reference_in = osr.SpatialReference()  
        pSpatial_reference_in.ImportFromEPSG(4326)    # WGS84 lat/lon
    else:
        pass

    if aAttribute_data is not None:
        aAttribute= aAttribute_data
        iFlag_attribute =1
    else:
        iFlag_attribute=0

    nVertex = len(aVertex_in)
    pDriver = ogr.GetDriverByName('GeoJSON')        
    pDataset_json = pDriver.CreateDataSource(sFilename_json_in)
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
        if iFlag_projected_in ==1:
            dummy1= Point( pVertex.dx, pVertex.dy )
            pass
        else:
            dummy1= Point( pVertex.dLongitude_degree, pVertex.dLatitude_degree ) 
            pass

        pGeometry_out = ogr.CreateGeometryFromWkb(dummy1.wkb)
        pFeature_out.SetGeometry(pGeometry_out)   
        pFeature_out.SetField("id", lID)
        if iFlag_attribute ==1:
            pFeature_out.SetField("con", int(aAttribute[i]) )
                
        pLayer_json.CreateFeature(pFeature_out)        
        lID =  lID + 1
        pass
        
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
        if iFlag_projected_in ==1:
            dummy1= Point( pVertex.dx, pVertex.dy )             
            pass
        else:
            dummy1= Point( pVertex.dLongitude_degree, pVertex.dLatitude_degree ) 
            pass

        pGeometry_out = ogr.CreateGeometryFromWkb(dummy1.wkb)
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


    


    
