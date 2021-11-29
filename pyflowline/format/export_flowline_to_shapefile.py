import os
import json
from osgeo import ogr, osr, gdal, gdalconst
from shapely.geometry import Point, LineString, MultiLineString
from shapely.wkt import loads
def export_flowline_to_shapefile(iFlag_projected_in, aFlowline_in, pSpatial_reference_in, \
    sFilename_shapefile_out, \
    aAttribute_field=None,\
    aAttribute_data=None,\
        aAttribute_dtype=None):
    """
    convert a shpefile to json format.
    This function should be used for stream flowline only.
    """

    if os.path.exists(sFilename_shapefile_out): 
        #delete it if it exists
        os.remove(sFilename_shapefile_out)
        pass

    nFlowline = len(aFlowline_in)

    if aAttribute_field is not None and aAttribute_data is not None and aAttribute_dtype is not None:
        iFlag_attribute = 1

        nAttribute1 = len(aAttribute_field)
        nAttribute2 = len(aAttribute_data)
        nAttribute3 = len(aAttribute_dtype)
        nAttribute4 = len(aAttribute_data[0])
        if nAttribute3 != nAttribute1 or nAttribute1 != nAttribute2 or nFlowline!= nAttribute4:
            print('The attribute is not correct, please check!')
            return
        else:
            iFlag_attribute = 1
    else:
        iFlag_attribute=0
        pass

    

    #pDriver_json = ogr.GetDriverByName('GeoJSON')
    pDriver_shapefile = ogr.GetDriverByName('ESRI Shapefile')
    #geojson
    pDataset_shapefile = pDriver_shapefile.CreateDataSource(sFilename_shapefile_out)  
    

    pLayer_shapefile = pDataset_shapefile.CreateLayer('flowline', pSpatial_reference_in, ogr.wkbLineString)
    # Add one attribute
    pLayer_shapefile.CreateField(ogr.FieldDefn('id', ogr.OFTInteger64)) #long type for high resolution

    #add the other fields
    if iFlag_attribute ==1:
        for i in range(nAttribute1):
            sField = aAttribute_field[i]
            dtype = aAttribute_dtype[i]
            if dtype == 'int':
                pLayer_shapefile.CreateField(ogr.FieldDefn(sField, ogr.OFTInteger64))
                pass
            else:
                pLayer_shapefile.CreateField(ogr.FieldDefn(sField, ogr.OFTReal))
                pass
        
    
    pLayerDefn = pLayer_shapefile.GetLayerDefn()
    pFeature_out = ogr.Feature(pLayerDefn)

    lID = 0
    for i in range(nFlowline):
        pFlowline = aFlowline_in[i]
        dummy =pFlowline.aVertex
        aPoint=list()
        for j in dummy:
            if iFlag_projected_in ==1:
                aPoint.append( Point( j.dx, j.dy ) )                
                pass
            else:
                aPoint.append( Point( j.dLongitude_degree, j.dLatitude_degree ) )
                pass

        dummy1= LineString( aPoint )
        pGeometry_out = ogr.CreateGeometryFromWkb(dummy1.wkb)
        pFeature_out.SetGeometry(pGeometry_out)
   
        pFeature_out.SetField("id", lID)
        if iFlag_attribute ==1:
            for k in range(nAttribute1):
                sField = aAttribute_field[k]
                dtype = aAttribute_dtype[k]
                dummy = aAttribute_data[k]
                if dtype == 'int':
                    pFeature_out.SetField(sField, int(dummy[i]))
                else:
                    pFeature_out.SetField(sField, float(dummy[i]))
        
        # Add new pFeature_shapefile to output Layer
        pLayer_shapefile.CreateFeature(pFeature_out)        
        lID =  lID + 1
        pass
        
    pDataset_shapefile.FlushCache()
    pDataset_shapefile = pLayer_shapefile = pFeature_out  = None    

    return


    
