
import os
import numpy as np
from osgeo import ogr, osr
from shapely.wkt import loads

from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_cell
from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_flowline

def intersect_flowline_with_mesh(iMesh_type_in, sFilename_mesh_in, sFilename_flowline_in, sFilename_output_in):

    if  os.path.exists(sFilename_mesh_in) and  os.path.exists(sFilename_flowline_in) : 
        pass
    else:
        print('The input file does not exist')
        return

    if os.path.exists(sFilename_output_in): 
        os.remove(sFilename_output_in)

    
    pDriver_geojson = ogr.GetDriverByName( "GeoJSON")
    aCell=list()
    aCell_intersect=list()
   
    pDataset_mesh = pDriver_geojson.Open(sFilename_mesh_in, 0)
    pDataset_flowline = pDriver_geojson.Open(sFilename_flowline_in, 0)   

    pLayer_mesh = pDataset_mesh.GetLayer(0)
    pSpatial_reference_mesh = pLayer_mesh.GetSpatialRef()
    nfeature_mesh = pLayer_mesh.GetFeatureCount()

    pLayer_flowline = pDataset_flowline.GetLayer(0)
    pSpatial_reference_flowline = pLayer_flowline.GetSpatialRef()
    nfeature_flowline = pLayer_flowline.GetFeatureCount()
    pLayerDefinition = pLayer_flowline.GetLayerDefn()
    
    comparison = pSpatial_reference_mesh.IsSame(pSpatial_reference_flowline)
    if(comparison != 1):
        iFlag_transform = 1
        transform = osr.CoordinateTransformation(pSpatial_reference_mesh, pSpatial_reference_flowline)
    else:
        iFlag_transform = 0

    pDataset_out = pDriver_geojson.CreateDataSource(sFilename_output_in)

    pLayerOut = pDataset_out.CreateLayer('flowline', pSpatial_reference_flowline, ogr.wkbMultiLineString)
    # Add one attribute
    pLayerOut.CreateField(ogr.FieldDefn('id', ogr.OFTInteger64)) #long type for high resolution
    pLayerOut.CreateField(ogr.FieldDefn('iseg', ogr.OFTInteger)) #long type for high resolution
    pLayerOut.CreateField(ogr.FieldDefn('iord', ogr.OFTInteger)) #long type for high resolution
    pLayerDefn = pLayerOut.GetLayerDefn()
    pFeatureOut = ogr.Feature(pLayerDefn)       
    lID_flowline = 0          
    aFlowline_intersect_all=list()   
    for pFeature_mesh in pLayer_mesh:       
        pGeometry_mesh = pFeature_mesh.GetGeometryRef()        
        dummy0 = loads( pGeometry_mesh.ExportToWkt() )
        aCoords_gcs = dummy0.exterior.coords
        aCoords_gcs= np.array(aCoords_gcs)       

        lCellID = pFeature_mesh.GetField("id")
        dLon = pFeature_mesh.GetField("lon")
        dLat = pFeature_mesh.GetField("lat")        
        dArea = pFeature_mesh.GetField("area")
        if (iFlag_transform ==1): 
            pGeometry_mesh.Transform(transform)
        if (pGeometry_mesh.IsValid()):
            pass
        else:
            print('Geometry issue')

        pGeometrytype_mesh = pGeometry_mesh.GetGeometryName()
        if(pGeometrytype_mesh == 'POLYGON'):            
            pCell = convert_gcs_coordinates_to_cell(iMesh_type_in, dLon, dLat, aCoords_gcs)     
            pCell.lCellID = lCellID             
            pCell.dArea = dArea 
            pCell.dLength = pCell.calculate_edge_length()
            pCell.dLength_flowline = pCell.dLength                     
            aFlowline_intersect = list()
            iFlag_intersected = 0 
            for j in range (nfeature_flowline):
                pFeature_flowline = pLayer_flowline.GetFeature(j)
                pGeometry_flowline = pFeature_flowline.GetGeometryRef()
                iStream_segment = pFeature_flowline.GetField("iseg")
                iStream_order = pFeature_flowline.GetField("iord")
                if (pGeometry_flowline.IsValid()):
                    pass
                else:
                    print('Geometry issue')         

                iFlag_intersect = pGeometry_flowline.Intersects( pGeometry_mesh )
                if( iFlag_intersect == True):
                    iFlag_intersected = 1
                    pGeometry_intersect = pGeometry_flowline.Intersection(pGeometry_mesh)                     
                    pGeometrytype_intersect = pGeometry_intersect.GetGeometryName()
                    if pGeometrytype_intersect == 'LINESTRING':
                        pFeatureOut.SetGeometry(pGeometry_intersect)
                        pFeatureOut.SetField("id", lID_flowline)         
                        pFeatureOut.SetField("iseg", iStream_segment)    
                        pFeatureOut.SetField("iord", iStream_order)           
                        pLayerOut.CreateFeature(pFeatureOut)    
                        dummy = loads( pGeometry_intersect.ExportToWkt() )
                        aCoords = dummy.coords                
                        dummy1= np.array(aCoords)
                        pLine = convert_gcs_coordinates_to_flowline(dummy1)
                        pLine.calculate_length()
                        pLine.lIndex = lID_flowline
                        pLine.iStream_segment = iStream_segment
                        pLine.iStream_order = iStream_order
                        aFlowline_intersect.append(pLine)
                        aFlowline_intersect_all.append(pLine)
                        lID_flowline = lID_flowline + 1
                    
                    else:
                        if(pGeometrytype_intersect == 'MULTILINESTRING'):
                            aLine = ogr.ForceToLineString(pGeometry_intersect)
                            for Line in aLine: 
                                pFeatureOut.SetGeometry(Line)
                                pFeatureOut.SetField("id", lID_flowline)         
                                pFeatureOut.SetField("iseg", iStream_segment)    
                                pFeatureOut.SetField("iord", iStream_order)           
                                pLayerOut.CreateFeature(pFeatureOut)    
                                dummy = loads( Line.ExportToWkt() )
                                aCoords = dummy.coords
                                dummy1= np.array(aCoords)
                                pLine = convert_gcs_coordinates_to_flowline(dummy1)
                                pLine.calculate_length()
                                pLine.lIndex = lID_flowline
                                pLine.iStream_segment = iStream_segment
                                pLine.iStream_order = iStream_order
                                aFlowline_intersect.append(pLine)
                                aFlowline_intersect_all.append(pLine)
                                lID_flowline = lID_flowline + 1
                            pass
                        else:
                            pass                            

                else:
                    pass
     
            #now add back to the cell object
            pCell.aFlowline = aFlowline_intersect
            pCell.nFlowline = len(aFlowline_intersect)
            if iFlag_intersected ==1:     
                pCell.iFlag_intersected = 1 
                pCell.dLength_flowline = 0.0
                for i in range (pCell.nFlowline):
                    pFlowline = pCell.aFlowline[i]
                    dLength_flowline = pFlowline.dLength 
                    if ( dLength_flowline > pCell.dLength_flowline ):
                        pCell.dLength_flowline = dLength_flowline

                #replace flowline length if there is an actual flowline                    
                aCell_intersect.append(pCell)
                aCell.append(pCell)
            else:
                pCell.iFlag_intersected = 0   
                aCell.append(pCell)
                pass
            
        else:
            pass

    
    return  aCell,aCell_intersect, aFlowline_intersect_all