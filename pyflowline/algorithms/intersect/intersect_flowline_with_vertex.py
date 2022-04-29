

import os
import numpy as np
from osgeo import ogr, osr
from shapely.wkt import loads

def intersect_flowline_with_vertex( sFilename_flowline_in, sFilename_vertex_in, sFilename_output_in):

    if  os.path.exists(sFilename_flowline_in) and  os.path.exists(sFilename_vertex_in) : 
        pass
    else:
        print('The input file does not exist')
        return

    if os.path.exists(sFilename_output_in): 

        os.remove(sFilename_output_in)

    
    pDriver_geojson = ogr.GetDriverByName( "GeoJSON")

    pDataset_flowline = pDriver_geojson.Open(sFilename_flowline_in, 0)
    pDataset_vertex = pDriver_geojson.Open(sFilename_vertex_in, 0)   

    pLayer_flowline = pDataset_flowline.GetLayer(0)
    pSpatial_reference_a = pLayer_flowline.GetSpatialRef()
    nfeature_flowline = pLayer_flowline.GetFeatureCount()

    pLayer_vertex = pDataset_vertex.GetLayer(0)
    pSpatial_reference_b = pLayer_vertex.GetSpatialRef()
    nfeature_vertex = pLayer_vertex.GetFeatureCount()
    pLayerDefinition = pLayer_vertex.GetLayerDefn()

    comparison = pSpatial_reference_a.IsSame(pSpatial_reference_b)
    if(comparison != 1):
        iFlag_transform = 1
        transform = osr.CoordinateTransformation(pSpatial_reference_a, pSpatial_reference_b)
    else:
        iFlag_transform = 0

    pDataset_out = pDriver_geojson.CreateDataSource(sFilename_output_in)

    pLayerOut = pDataset_out.CreateLayer('flowline', pSpatial_reference_b, ogr.wkbMultiPoint)
    # Add one attribute
    pLayerOut.CreateField(ogr.FieldDefn('id', ogr.OFTInteger64)) #long type for high resolution
    
    pLayerDefn = pLayerOut.GetLayerDefn()
    pFeatureOut = ogr.Feature(pLayerDefn)    
   
    lID_flowline = 0           

    aVertex_intersect=list()
    #for i in range (nfeature_mesh):
    for pFeature_flowline in pLayer_flowline:
       
        #pFeature_mesh= pLayer_mesh.GetFeature(i)
        pGeometry_flowline = pFeature_flowline.GetGeometryRef()        
        dummy0 = loads( pGeometry_flowline.ExportToWkt() )
        aCoords_gcs = dummy0.coords
        aCoords_gcs= np.array(aCoords_gcs)       

        if (iFlag_transform ==1): #projections are different
            pGeometry_flowline.Transform(transform)

        if (pGeometry_flowline.IsValid()):
            pass
        else:
            print('Geometry issue')

        #convert geometry to edge
        pGeometrytype_flowline = pGeometry_flowline.GetGeometryName()
        if(pGeometrytype_flowline == 'LINESTRING'):               
                     
            aFlowline_intersect = list()
            iFlag_intersected = 0 
            for j in range (nfeature_vertex):
            #for pFeature_flowline in pLayer_flowline:
                pFeature_vertex = pLayer_vertex.GetFeature(j)
                pGeometry_vertex = pFeature_vertex.GetGeometryRef()
                if (pGeometry_vertex.IsValid()):
                    pass
                else:
                    print('Geometry issue')

                iFlag_intersect = pGeometry_vertex.Intersects( pGeometry_flowline )
                if( iFlag_intersect == True):
                    iFlag_intersected = 1
                    pGeometry_intersect = pGeometry_vertex.Intersection(pGeometry_flowline)      
                    #add more process here to 
                    pGeometrytype_intersect = pGeometry_intersect.GetGeometryName()
                    if pGeometrytype_intersect == 'LINESTRING':
                        npoint = pGeometry_intersect.GetGeometryCount()
                        for i  in range(npoint): 
                            point = pGeometry_intersect.GetGeometryRef(i) 

                    else:
                        pass
                             
                    
                else:
                    pass
           
            
        else:
            pass

    pDataset_out = pLayerOut = None    

    return   aVertex_intersect