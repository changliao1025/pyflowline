
import numpy as np
from osgeo import ogr, osr, gdal
from shapely.wkt import loads

from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_flowline

def read_nhdplus_flowline_shapefile_attribute(sFilename_shapefile_in):
    """
    convert a shpefile to json format.
    This function should be used for stream flowline only.
    """

    #aFromNode=list()
    #aToNode=list()
    aNHDPlusID=list()

    #pDriver_json = ogr.GetDriverByName('GeoJSON')
    pDriver_shapefile = ogr.GetDriverByName('ESRI Shapefile')
   
    pDataset_shapefile = pDriver_shapefile.Open(sFilename_shapefile_in, gdal.GA_ReadOnly)
    pLayer_shapefile = pDataset_shapefile.GetLayer(0)
    pSpatialRef_shapefile = pLayer_shapefile.GetSpatialRef()    
        
    for pFeature_shapefile in pLayer_shapefile:
        
        pGeometry_in = pFeature_shapefile.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()

        #sFromNode = str(pFeature_shapefile.GetField("NHDPlusF_4"))
        #sToNode = str(pFeature_shapefile.GetField("NHDPlusF_5"))
        #if (sFromNode !='None')  & ( sToNode != 'None' ):
            #lFromNode = int(float(sFromNode))
            #lToNode = int(float(sToNode))
            #lFromNode = int(pFeature_shapefile.GetField("NHDPlusF_4"))
            #lToNode = int(pFeature_shapefile.GetField("NHDPlusF_4"))
        lNHDPlusID = int(pFeature_shapefile.GetField("NHDPlusID"))
            #aFromNode.append(lFromNode)
            #aToNode.append(lToNode)
        aNHDPlusID.append(lNHDPlusID)
        
    
    #we also need to spatial reference
    #aFromNode, aToNode, 
    return aNHDPlusID

def extract_nhdplus_flowline_shapefile_by_attribute(sFilename_shapefile_in, aAttribute_in):
    aFlowline =list()
    #pDriver_json = ogr.GetDriverByName('GeoJSON')
    pDriver_shapefile = ogr.GetDriverByName('ESRI Shapefile')
   
    pDataset_shapefile = pDriver_shapefile.Open(sFilename_shapefile_in, gdal.GA_ReadOnly)
    pLayer_shapefile = pDataset_shapefile.GetLayer(0)
    pSpatialRef_shapefile = pLayer_shapefile.GetSpatialRef()

    pSpatial_reference_gcs = osr.SpatialReference()
    pSpatial_reference_gcs.ImportFromEPSG(4326)
    pSpatial_reference_gcs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

    comparison = pSpatialRef_shapefile.IsSame(pSpatial_reference_gcs)
    if(comparison != 1):
        iFlag_transform =1
        pTransform = osr.CoordinateTransformation(pSpatialRef_shapefile, pSpatial_reference_gcs)
    else:
        iFlag_transform =0   

    lID = 0
    for pFeature_shapefile in pLayer_shapefile:
        
        pGeometry_in = pFeature_shapefile.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()

        #sFromNode = str(pFeature_shapefile.GetField("NHDPlusF_4"))
        #sToNode = str(pFeature_shapefile.GetField("NHDPlusF_5"))
        #if (sFromNode !='None')  & ( sToNode != 'None' ):
            #lFromNode = int(float(sFromNode))
            #lToNode = int(float(sToNode))
            #lFromNode = int(pFeature_shapefile.GetField("NHDPlusF_4"))
            #lToNode = int(pFeature_shapefile.GetField("NHDPlusF_4"))
        lNHDPlusID = int(pFeature_shapefile.GetField("NHDPlusID"))
            
        if (iFlag_transform ==1): #projections are different
            pGeometry_in.Transform(pTransform)
        if (pGeometry_in.IsValid()):
            pass
        else:
            print('Geometry issue')
        if lNHDPlusID in aAttribute_in:
            if(sGeometry_type == 'MULTILINESTRING'):
                aLine = ogr.ForceToLineString(pGeometry_in)
                for Line in aLine: 
                    dummy = loads( Line.ExportToWkt() )
                    aCoords = dummy.coords
                    #pLine= LineString( aCoords[::-1 ] )
                    dummy1= np.array(aCoords)
                    pLine = convert_gcs_coordinates_to_flowline(dummy1)
                    pLine.lIndex = lID                    
                    aFlowline.append(pLine)
                    lID = lID + 1
            else:
                if sGeometry_type =='LINESTRING':
                    dummy = loads( pGeometry_in.ExportToWkt() )
                    aCoords = dummy.coords
                    #pLine= LineString( aCoords[::-1 ] )
                    dummy1= np.array(aCoords)
                    pLine = convert_gcs_coordinates_to_flowline(dummy1)
                    pLine.lIndex = lID                    
                    aFlowline.append(pLine)
                    lID = lID + 1
                else:
                    print(sGeometry_type)
                    pass
        else:
            pass

    return aFlowline

def track_nhdplus_flowline(aNHDPlusID_filter_in, aFromFlowline_in, aToFlowline_in, lNHDPlusID_in):
    #aNHDPlusID_dam_headwater = list()
    aNHDPlusID_dam_nonheadwater = list()
    def tag_downstream(lNHDPlusID_from):
        if lNHDPlusID_from in aNHDPlusID_filter_in:
            dummy_index = aNHDPlusID_filter_in.index(lNHDPlusID_from)
            pass
        else:
            if lNHDPlusID_from in aFromFlowline_in:

                if lNHDPlusID_from in aNHDPlusID_dam_nonheadwater:
                    pass
                else:
                    aNHDPlusID_dam_nonheadwater.append(lNHDPlusID_from)
                    dummy_index = np.where(aFromFlowline_in == lNHDPlusID_from )
                    nDownstream = dummy_index[0].size
                    for i in range(nDownstream):
                        lNHDPlusID_to = aToFlowline_in[dummy_index[0][i]   ]                    
                        if lNHDPlusID_to==0:
                            pass
                        else:
                            #print(lNHDPlusID_in, lNHDPlusID_to)
                            tag_downstream(lNHDPlusID_to)
        return 

    if lNHDPlusID_in in aNHDPlusID_filter_in:
        dummy_index = aNHDPlusID_filter_in.index(lNHDPlusID_in)
        #print( lNHDPlusID, dummy_index)        
    else:
        tag_downstream(lNHDPlusID_in)   
        #remove the first one
        aNHDPlusID_dam_nonheadwater.pop(0)    


    return  aNHDPlusID_dam_nonheadwater