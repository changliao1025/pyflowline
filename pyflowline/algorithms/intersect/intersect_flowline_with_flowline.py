import os
import numpy as np
from osgeo import ogr, osr
from shapely.wkt import loads
from pyflowline.classes.vertex import pyvertex
from pyflowline.algorithms.auxiliary.find_index_in_list import find_vertex_in_list

def intersect_flowline_with_flowline( sFilename_flowline_a_in, sFilename_flowline_b_in, sFilename_output_in):
    if  os.path.exists(sFilename_flowline_a_in) and  os.path.exists(sFilename_flowline_b_in) : 
        pass
    else:
        print('The input file does not exist')
        return

    if os.path.exists(sFilename_output_in): 
        os.remove(sFilename_output_in)
    
    pDriver_geojson = ogr.GetDriverByName( "GeoJSON")
    pDataset_flowline_a = pDriver_geojson.Open(sFilename_flowline_a_in, 0)
    pDataset_flowline_b = pDriver_geojson.Open(sFilename_flowline_b_in, 0)   
    pLayer_flowline_a = pDataset_flowline_a.GetLayer(0)
    pSpatial_reference_a = pLayer_flowline_a.GetSpatialRef()
    nfeature_flowline_a = pLayer_flowline_a.GetFeatureCount()
    pLayer_flowline_b = pDataset_flowline_b.GetLayer(0)
    pSpatial_reference_b = pLayer_flowline_b.GetSpatialRef()
    nfeature_flowline_b = pLayer_flowline_b.GetFeatureCount()
    pLayerDefinition = pLayer_flowline_b.GetLayerDefn()
    schema =list()
    for n in range(pLayerDefinition.GetFieldCount()):
        fdefn = pLayerDefinition.GetFieldDefn(n)
        schema.append(fdefn.name)
    if 'id' in schema:
        iFlag_id = 1
    else:
        iFlag_id = 0   
    
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
   
    lVertexID = 0           

    aVertex_intersect=list()
    #for i in range (nfeature_mesh):
    for pFeature_flowline_a in pLayer_flowline_a:
       
        #pFeature_mesh= pLayer_mesh.GetFeature(i)
        pGeometry_flowline_a = pFeature_flowline_a.GetGeometryRef()        
        dummy0 = loads( pGeometry_flowline_a.ExportToWkt() )
        aCoords_gcs = dummy0.coords
        aCoords_gcs= np.array(aCoords_gcs)       

        if (iFlag_transform ==1): #projections are different
            pGeometry_flowline_a.Transform(transform)

        if (pGeometry_flowline_a.IsValid()):
            pass
        else:
            print('Geometry issue')

        #convert geometry to edge
        pGeometrytype_flowline_a = pGeometry_flowline_a.GetGeometryName()
        if(pGeometrytype_flowline_a == 'LINESTRING'):            
            
                     
            aFlowline_intersect = list()
            iFlag_intersected = 0 
            for j in range (nfeature_flowline_b):
            #for pFeature_flowline in pLayer_flowline:
                pFeature_flowline_b = pLayer_flowline_b.GetFeature(j)
                pGeometry_flowline_b = pFeature_flowline_b.GetGeometryRef()

                if iFlag_id ==1:
                    lFlowlineID = pFeature_flowline_b.GetField("id")
                else:
                    lFlowlineID = -1

                if (pGeometry_flowline_b.IsValid()):
                    pass
                else:
                    print('Geometry issue')

                iFlag_intersect = pGeometry_flowline_b.Intersects( pGeometry_flowline_a )
                if( iFlag_intersect == True):
                    iFlag_intersected = 1
                    pGeometry_intersect = pGeometry_flowline_b.Intersection(pGeometry_flowline_a)      
                    #add more process here to 
                    pGeometrytype_intersect = pGeometry_intersect.GetGeometryName()
                    
                    if pGeometrytype_intersect == 'MULTIPOINT':
                        npoint = pGeometry_intersect.GetGeometryCount()
                        for i  in range(npoint): 
                            point = pGeometry_intersect.GetGeometryRef(i)                            

                            point0= dict()   
                            point0['dLongitude_degree'] = point.GetX()
                            point0['dLatitude_degree'] = point.GetY()
                            pVertex=pyvertex(point0)
                            pVertex.lFlowlineID = lFlowlineID
                            iFlag_exist, lIndex = find_vertex_in_list( aVertex_intersect,  pVertex)
                            if iFlag_exist ==1:
                                pass
                            else:                       
                                aVertex_intersect.append(pVertex)
                                pFeatureOut.SetGeometry(point)
                                pFeatureOut.SetField("id", lVertexID)         
                                pLayerOut.CreateFeature(pFeatureOut)    
                                lVertexID = lVertexID + 1

                    else:
                        #this branch possible has error, report an issue if crash
                        point= dict()   
                        point['dLongitude_degree'] = pGeometry_intersect.GetX()
                        point['dLatitude_degree'] = pGeometry_intersect.GetY()
                        pVertex=pyvertex(point)
                        pVertex.lFlowlineID = lFlowlineID
                        iFlag_exist, lIndex = find_vertex_in_list( aVertex_intersect,  pVertex)
                        if iFlag_exist ==1:
                            pass
                        else:                
                            aVertex_intersect.append(pVertex)
                            pFeatureOut.SetGeometry(pGeometry_intersect)
                            pFeatureOut.SetField("id", lVertexID)         
                            pLayerOut.CreateFeature(pFeatureOut)    
                            lVertexID = lVertexID + 1                                        
                    
                else:
                    pass           
            
        else:
            pass

    pDataset_out = pLayerOut = None    

    return   aVertex_intersect