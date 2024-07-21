
import os
import importlib.util
import numpy as np
from osgeo import ogr, osr
#from shapely.wkt import loads

from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_cell
from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_flowline
from pyearth.gis.location.get_geometry_coordinates import get_geometry_coordinates

iFlag_cython = importlib.util.find_spec("cython")
if iFlag_cython is not None:
    from tinyr import RTree
    iFlag_use_rtree = 1
    from pyflowline.algorithms.cython.kernel import find_vertex_in_list
else:
    iFlag_use_rtree = 0
    from pyflowline.algorithms.auxiliary.find_vertex_in_list import find_vertex_in_list
    pass

def intersect_flowline_with_mesh(iMesh_type_in, sFilename_mesh_in, sFilename_flowline_in, sFilename_output_in):
    if  os.path.exists(sFilename_mesh_in) and  os.path.exists(sFilename_flowline_in) :
        pass
    else:
        print('The input file does not exist')
        return

    if os.path.exists(sFilename_output_in):
        os.remove(sFilename_output_in)


    pDriver_geojson = ogr.GetDriverByName( "GeoJSON")
    #aCell=list()
    aCell_intersect=list()

    pDataset_mesh = pDriver_geojson.Open(sFilename_mesh_in, 0)
    pLayer_mesh = pDataset_mesh.GetLayer(0)
    pSpatial_reference_mesh = pLayer_mesh.GetSpatialRef()
    nfeature_mesh = pLayer_mesh.GetFeatureCount()

    pDataset_flowline = pDriver_geojson.Open(sFilename_flowline_in, 0)
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
    pLayerOut.CreateField(ogr.FieldDefn('lineid', ogr.OFTInteger64)) #long type for high resolution
    pLayerOut.CreateField(ogr.FieldDefn('stream_segment', ogr.OFTInteger)) #long type for high resolution
    pLayerOut.CreateField(ogr.FieldDefn('stream_order', ogr.OFTInteger)) #long type for high resolution
    pLayerDefn = pLayerOut.GetLayerDefn()
    pFeatureOut = ogr.Feature(pLayerDefn)

    lFlowlineID = 1
    aFlowline_intersect_all=list()
    if iFlag_use_rtree ==1: #use the rtree to speed up
        #index_flowline = rtree.index.Index()
        interleaved = True
        index_flowline = RTree(interleaved=interleaved, max_cap=5, min_cap=2)
        for i in range(nfeature_flowline):
            lID = i
            pFeature_flowline = pLayer_flowline.GetFeature(i)
            pGeometry_flowline = pFeature_flowline.GetGeometryRef()
            left, right, bottom, top= pGeometry_flowline.GetEnvelope()
            pBound= (left, bottom, right, top)
            index_flowline.insert(lID, pBound)  #

        #now intersect using rtree
        for j in range (nfeature_mesh):
            pFeature_mesh = pLayer_mesh.GetFeature(j)
            pGeometry_mesh = pFeature_mesh.GetGeometryRef()
            aCoords_gcs = get_geometry_coordinates(pGeometry_mesh)
            lCellID = pFeature_mesh.GetField("cellid")
            dLon = pFeature_mesh.GetField("longitude")
            dLat = pFeature_mesh.GetField("latitude")
            dArea = pFeature_mesh.GetField("area")
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

                left, right, bottom, top= pGeometry_mesh.GetEnvelope()
                pBound= (left, bottom, right, top)

                aIntersect = list(index_flowline.search(pBound))
                for k in aIntersect:
                    pFeature_flowline = pLayer_flowline.GetFeature(k)
                    pGeometry_flowline = pFeature_flowline.GetGeometryRef()
                    iFlag_intersect = pGeometry_flowline.Intersects( pGeometry_mesh )
                    if( iFlag_intersect == True):
                        iFlag_intersected = 1
                        pGeometry_intersect = pGeometry_flowline.Intersection(pGeometry_mesh)
                        pGeometrytype_intersect = pGeometry_intersect.GetGeometryName()
                        iStream_segment = pFeature_flowline.GetField("stream_segment")
                        iStream_order = pFeature_flowline.GetField("stream_order")
                        if pGeometrytype_intersect == 'LINESTRING':
                            pFeatureOut.SetGeometry(pGeometry_intersect)
                            pFeatureOut.SetField("lineid", lFlowlineID)
                            pFeatureOut.SetField("stream_segment", iStream_segment)
                            pFeatureOut.SetField("stream_order", iStream_order)
                            pLayerOut.CreateFeature(pFeatureOut)

                            aCoords = list()
                            for i in range(0, pGeometry_intersect.GetPointCount()):
                                pt = pGeometry_intersect.GetPoint(i)
                                aCoords.append( [ pt[0], pt[1]])

                            dummy1= np.array(aCoords)
                            pFlowline = convert_gcs_coordinates_to_flowline(dummy1)
                            pFlowline.calculate_length()
                            pFlowline.lFlowlineIndex = lFlowlineID
                            pFlowline.iStream_segment = iStream_segment
                            pFlowline.iStream_order = iStream_order
                            aFlowline_intersect.append(pFlowline)
                            aFlowline_intersect_all.append(pFlowline)
                            lFlowlineID = lFlowlineID + 1

                        else:
                            if(pGeometrytype_intersect == 'MULTILINESTRING'):
                                nLine = pGeometry_intersect.GetGeometryCount()
                                for i in range(nLine):
                                    Line = pGeometry_intersect.GetGeometryRef(i)
                                    pFeatureOut.SetGeometry(Line)
                                    pFeatureOut.SetField("lineid", lFlowlineID)
                                    pFeatureOut.SetField("stream_segment", iStream_segment)
                                    pFeatureOut.SetField("stream_order", iStream_order)
                                    pLayerOut.CreateFeature(pFeatureOut)
                                    aCoords = list()
                                    for i in range(0, Line.GetPointCount()):
                                        pt = Line.GetPoint(i)
                                        aCoords.append( [ pt[0], pt[1]])

                                    dummy1= np.array(aCoords)
                                    pFlowline = convert_gcs_coordinates_to_flowline(dummy1)
                                    pFlowline.calculate_length()
                                    pFlowline.lFlowlineIndex = lFlowlineID
                                    pFlowline.iStream_segment = iStream_segment
                                    pFlowline.iStream_order = iStream_order
                                    aFlowline_intersect.append(pFlowline)
                                    aFlowline_intersect_all.append(pFlowline)
                                    lFlowlineID = lFlowlineID + 1
                                pass
                            else:
                                pass


                        pass

                #now add back to the cell object
                pCell.aFlowline = aFlowline_intersect
                pCell.nFlowline = len(aFlowline_intersect)
                if iFlag_intersected ==1:
                    pCell.iFlag_intersected = 1
                    pCell.dLength_flowline = 0.0 #reset the flowline length
                    for i in range (pCell.nFlowline):
                        pFlowline = pCell.aFlowline[i]
                        dLength_flowline = pFlowline.dLength
                        if ( dLength_flowline > pCell.dLength_flowline ):
                            pCell.dLength_flowline = dLength_flowline

                    #replace flowline length if there is an actual flowline
                    aCell_intersect.append(pCell)
                    #aCell.append(pCell)
                else:
                    pCell.iFlag_intersected = 0
                    #aCell.append(pCell)
                    pass
    else:

        for pFeature_mesh in pLayer_mesh:
            pGeometry_mesh = pFeature_mesh.GetGeometryRef()
            aCoords_gcs = get_geometry_coordinates(pGeometry_mesh)

            lCellID = pFeature_mesh.GetField("cellid")
            dLon = pFeature_mesh.GetField("longitude")
            dLat = pFeature_mesh.GetField("latitude")
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
                    iStream_segment = pFeature_flowline.GetField("stream_segment")
                    iStream_order = pFeature_flowline.GetField("stream_order")
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
                            pFeatureOut.SetField("lineid", lFlowlineID)
                            pFeatureOut.SetField("stream_segment", iStream_segment)
                            pFeatureOut.SetField("stream_order", iStream_order)
                            pLayerOut.CreateFeature(pFeatureOut)

                            aCoords = list()
                            for i in range(0, pGeometry_intersect.GetPointCount()):
                                pt = pGeometry_intersect.GetPoint(i)
                                aCoords.append( [ pt[0], pt[1]])

                            dummy1= np.array(aCoords)
                            pFlowline = convert_gcs_coordinates_to_flowline(dummy1)
                            pFlowline.calculate_length()
                            pFlowline.lFlowlineIndex = lFlowlineID
                            pFlowline.iStream_segment = iStream_segment
                            pFlowline.iStream_order = iStream_order
                            aFlowline_intersect.append(pFlowline)
                            aFlowline_intersect_all.append(pFlowline)
                            lFlowlineID = lFlowlineID + 1

                        else:
                            if(pGeometrytype_intersect == 'MULTILINESTRING'):
                                nLine = pGeometry_intersect.GetGeometryCount()
                                for i in range(nLine):
                                    Line = pGeometry_intersect.GetGeometryRef(i)
                                    pFeatureOut.SetGeometry(Line)
                                    pFeatureOut.SetField("lineid", lFlowlineID)
                                    pFeatureOut.SetField("stream_segment", iStream_segment)
                                    pFeatureOut.SetField("stream_order", iStream_order)
                                    pLayerOut.CreateFeature(pFeatureOut)
                                    aCoords = list()
                                    for i in range(0, Line.GetPointCount()):
                                        pt = Line.GetPoint(i)
                                        aCoords.append( [ pt[0], pt[1]])

                                    dummy1= np.array(aCoords)
                                    pFlowline = convert_gcs_coordinates_to_flowline(dummy1)
                                    pFlowline.calculate_length()
                                    pFlowline.lFlowlineIndex = lFlowlineID
                                    pFlowline.iStream_segment = iStream_segment
                                    pFlowline.iStream_order = iStream_order
                                    aFlowline_intersect.append(pFlowline)
                                    aFlowline_intersect_all.append(pFlowline)
                                    lFlowlineID = lFlowlineID + 1
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
                    #aCell.append(pCell)
                else:
                    pCell.iFlag_intersected = 0
                    #aCell.append(pCell)
                    pass

            else:
                pass



    #quality control
    #find the flowline that has the largest stream segment index
    pVertex_outlet = None
    iSegment_max = 0
    iIndex_flowline = -1
    for i in range(nfeature_flowline):
        pFeature_flowline = pLayer_flowline.GetFeature(i)
        iStream_segment = pFeature_flowline.GetField("stream_segment")
        if (iStream_segment > iSegment_max):
            iSegment_max = iStream_segment
            iIndex_flowline = i

    if (iIndex_flowline != -1):
        #find all the flowlines and cells that have the same stream segment

        #get the first vertex of the flowline
        pFeature_flowline = pLayer_flowline.GetFeature(iIndex_flowline)
        pGeometry_flowline = pFeature_flowline.GetGeometryRef()
        aCoords = list()
        for i in range(0, pGeometry_flowline.GetPointCount()):
            pt = pGeometry_flowline.GetPoint(i)
            aCoords.append( [ pt[0], pt[1]])
        dummy1= np.array(aCoords)
        pFlowline = convert_gcs_coordinates_to_flowline(dummy1)
        pVertex_start = pFlowline.pVertex_start
        pVertex_outlet = pFlowline.pVertex_end
        aFlowline_last = list()
        nFlowline_intersect = len(aFlowline_intersect_all)
        for i in range(nFlowline_intersect):
            pFlowline = aFlowline_intersect_all[i]
            if ( pFlowline.iStream_segment == iSegment_max):
                aFlowline_last.append(pFlowline)


        aFlowline_stay = list()
        #now check the connectivity of the flowlines because the flowlines are not necessarily connected
        iFlag_done = 0
        while (iFlag_done ==0):
            iFlag_found = 0
            nFlowline_last = len(aFlowline_last)
            for i in range(nFlowline_last):
                pFlowline = aFlowline_last[i]
                pVertex_start_dummy = pFlowline.pVertex_start
                pVertex_end_dummy = pFlowline.pVertex_end
                if pVertex_start_dummy == pVertex_start:
                    #found the head
                    iFlag_found = 1
                    pVertex_start = pVertex_end_dummy
                    aFlowline_stay.append(i)
                    pVertex_outlet = pVertex_end_dummy
                    break

            if iFlag_found == 1:
                iFlag_done = 0
            else:
                iFlag_done = 1

        #check the size of stay
        #nFlowline_stay = len(aFlowline_stay)
        #if (nFlowline_stay == nFlowline_intersect):
        #    print('All flowlines in the last segment are connected')
        #else:
        #    for i in range(nFlowline_last):
        #        if i not in aFlowline_stay:
        #            pFlowline = aFlowline_last[i]
        #            #remove this flowline from the result
        #            pass


    #return  aCell, aCell_intersect, aFlowline_intersect_all, pVertex_outlet
    return   aCell_intersect, aFlowline_intersect_all, pVertex_outlet