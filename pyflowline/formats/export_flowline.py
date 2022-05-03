import os
import json
from osgeo import ogr, osr
from shapely.geometry import Point, LineString
from pyflowline.classes.edge import pyedge
from pyflowline.classes.link import pycelllink

def export_flowline_to_json( aFlowline_in, \
    sFilename_json_in, \
    iFlag_projected_in= None, \
    pSpatial_reference_in=None, \
    aAttribute_field=None,\
    aAttribute_data=None,\
    aAttribute_dtype=None):
    """
    convert a shpefile to json format.
    This function should be used for stream flowline only.
    """

    if os.path.exists(sFilename_json_in): 
        os.remove(sFilename_json_in)
        pass

    nFlowline = len(aFlowline_in)
    if iFlag_projected_in is None:
        iFlag_projected_in = 0
    else:
        iFlag_projected_in = 1

    if  pSpatial_reference_in is None:        
        pSpatial_reference_in = osr.SpatialReference()  
        pSpatial_reference_in.ImportFromEPSG(4326)    # WGS84 lat/lon
    else:
        pass

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


    pDriver_json = ogr.GetDriverByName('GeoJSON')
    
    pDataset_json = pDriver_json.CreateDataSource(sFilename_json_in)  
    

    pLayer_json = pDataset_json.CreateLayer('flowline', pSpatial_reference_in, ogr.wkbLineString)
    # Add one attribute
    pLayer_json.CreateField(ogr.FieldDefn('id', ogr.OFTInteger64)) #long type for high resolution

    #add the other fields
    if iFlag_attribute ==1:
        for i in range(nAttribute1):
            sField = aAttribute_field[i]
            dtype = aAttribute_dtype[i]
            if dtype == 'int':
                pLayer_json.CreateField(ogr.FieldDefn(sField, ogr.OFTInteger64))
                pass
            else:
                pLayer_json.CreateField(ogr.FieldDefn(sField, ogr.OFTReal))
                pass
        
    
    pLayerDefn = pLayer_json.GetLayerDefn()
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
        pLayer_json.CreateFeature(pFeature_out)        
        lID =  lID + 1
        pass
        
    pDataset_json.FlushCache()
    pDataset_json = pLayer_json = pFeature_out  = None    

    return

    

def export_flowline_to_shapefile(iFlag_projected_in, aFlowline_in, pSpatial_reference_in, \
    sFilename_shapefile_in, \
    aAttribute_field=None,\
    aAttribute_data=None,\
        aAttribute_dtype=None):
    """
    convert a shpefile to json format.
    This function should be used for stream flowline only.
    """

    if os.path.exists(sFilename_shapefile_in): 
        #delete it if it exists
        os.remove(sFilename_shapefile_in)
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

    pDriver_shapefile = ogr.GetDriverByName('ESRI Shapefile')
    pDataset_shapefile = pDriver_shapefile.CreateDataSource(sFilename_shapefile_in)      
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

  

def export_flowline_info_to_json(aCell, aCell_intersect_in, aFlowline_in, sFilename_json_out):    
    #export the flowline topology to json
    ncell= len(aCell_intersect_in)
    nflowline= len(aFlowline_in)
    aLink =list()   
    for i in range(1, nflowline+1):
        pFlowline = aFlowline_in[i-1]
        nVertex = pFlowline.nVertex
        nEdge = pFlowline.nEdge
        for j in range(1, nEdge+1):         
            pEdge = pFlowline.aEdge[j-1]
            pVertex_start = pEdge.pVertex_start
            pVertex_end = pEdge.pVertex_end
            for k in range(ncell):
                if aCell_intersect_in[k].pVertex_center == pVertex_start:
                    pMpas_start = aCell_intersect_in[k]
                    pass
                if aCell_intersect_in[k].pVertex_center == pVertex_end:
                    pMpas_end = aCell_intersect_in[k]
                    pass
            pEdge_link = pyedge(pVertex_start, pVertex_end)  
            pLink = pycelllink(pMpas_start, pMpas_end, pEdge_link)
            aLink.append(pLink)

    with open(sFilename_json_out, 'w', encoding='utf-8') as f:
        sJson = json.dumps([json.loads(ob.tojson()) for ob in aLink], indent = 4)
        f.write(sJson)    
        f.close()
                
    
    return

    
   