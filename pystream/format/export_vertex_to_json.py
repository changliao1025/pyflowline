import os
import json
from osgeo import ogr, osr, gdal, gdalconst
from shapely.geometry import Point, LineString, MultiLineString
from shapely.wkt import loads
def export_vertex_to_json(aVertex_in, pSpatial_reference_in, sFilename_json_out):
    """
    convert a shpefile to json format.
    This function should be used for stream flowline only.
    """

    if os.path.exists(sFilename_json_out): 
        #delete it if it exists
        os.remove(sFilename_json_out)
        pass

    pDriver = ogr.GetDriverByName('GeoJSON')
    #geojson
    pDataset_json = pDriver.CreateDataSource(sFilename_json_out)

   
    

    pLayer_json = pDataset_json.CreateLayer('vertex', pSpatial_reference_in, ogr.wkbMultiLineString)
    # Add one attribute
    pLayer_json.CreateField(ogr.FieldDefn('id', ogr.OFTInteger64)) #long type for high resolution
    
    pLayerDefn = pLayer_json.GetLayerDefn()
    pFeature_out = ogr.Feature(pLayerDefn)

    lID = 0
    for pVertex in aVertex_in:       

        dummy1= Point( pVertex.dx, pVertex.dy )
        pGeometry_out = ogr.CreateGeometryFromWkb(dummy1.wkb)
        pFeature_out.SetGeometry(pGeometry_out)
   
        pFeature_out.SetField("id", lID)
        
        # Add new pFeature_shapefile to output Layer
        pLayer_json.CreateFeature(pFeature_out)        
        lID =  lID + 1
        pass
        
    pDataset_json.FlushCache()
    pDataset_json = pLayer_json = pFeature_out  = None    

    return


    
