import os
import numpy as np
from osgeo import ogr, osr
from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_cell
from pyearth.gis.spatialref.reproject_coordinates import reproject_coordinates_batch

def create_triangular_mesh(dX_left_in, dY_bot_in,
                    dResolution_meter_in,
                    ncolumn_in, nrow_in,
                    sFilename_output_in,
                    pProjection_reference_in,
                    pBoundary_in= None):

    #for the reason that a geometry object will be crash if the associated dataset is closed, we must pass wkt string
    #https://gdal.org/api/python_gotchas.html
    if pBoundary_in is None:
        print('Are you sure you want to create a mesh without boundary?')
        print('In this case, the program will generate a boundary for you using the bounding box of the mesh')
        #create a boundary for the mesh
        pBoundary = ogr.Geometry(ogr.wkbPolygon)
        ring = ogr.Geometry(ogr.wkbLinearRing)
        x1 = dX_left_in
        y1 = dY_bot_in
        x2 = x1 + ncolumn_in * dResolution_meter_in
        y2 = y1
        x3 = x2
        y3 = y1 + nrow_in * dResolution_meter_in
        x4 = x1
        y4 = y3
        ring.AddPoint(x1, y1)
        ring.AddPoint(x2, y2)
        ring.AddPoint(x3, y3)
        ring.AddPoint(x4, y4)
        ring.AddPoint(x1, y1)
        pBoundary.AddGeometry(ring)
    else:
        pBoundary = ogr.CreateGeometryFromWkt(pBoundary_in)

    if os.path.exists(sFilename_output_in):
        #delete it if it exists
        os.remove(sFilename_output_in)

    pDriver_geojson = ogr.GetDriverByName('GeoJSON')

    pSpatial_reference_gcs = osr.SpatialReference()
    pSpatial_reference_gcs.ImportFromEPSG(4326)    # WGS84 lat/lon

    pDataset = pDriver_geojson.CreateDataSource(sFilename_output_in)
    pLayer = pDataset.CreateLayer('cell', pSpatial_reference_gcs, ogr.wkbPolygon)
    # Add one attribute
    pLayer.CreateField(ogr.FieldDefn('cellid', ogr.OFTInteger64)) #long type for high resolution
    pLayer.CreateField(ogr.FieldDefn('longitude', ogr.OFTReal)) #long type for high resolution
    pLayer.CreateField(ogr.FieldDefn('latitude', ogr.OFTReal)) #long type for high resolution
    pArea_field = ogr.FieldDefn('area', ogr.OFTReal)
    pArea_field.SetWidth(20)
    pArea_field.SetPrecision(2)
    pLayer.CreateField(pArea_field)
    pLayerDefn = pLayer.GetLayerDefn()
    pFeature = ogr.Feature(pLayerDefn)
    xleft = dX_left_in
    ybottom = dY_bot_in
    dArea = np.power(dResolution_meter_in,2.0)
    #tin edge
    dLength_edge = np.sqrt(  4.0 * dArea /  np.sqrt(3.0) )
    dX_shift = 0.5 * dLength_edge
    dY_shift = 0.5 * dLength_edge * np.sqrt(3.0)
    dX_spacing = dX_shift * 2
    dY_spacing = dY_shift


    lCellID = 1

    #geojson
    aTriangular=list()
    #.........
    #(x2,y2)-----(x3,y3)
    #   |           |
    #(x1,y1)-----(x4,y4)
    #...............
    for column in range(0, ncolumn_in):
        for row in range(0, nrow_in):

            if column % 2 == 0 :
                if row % 2 == 0:
                    #define a polygon here
                    x1 = xleft + (column * dX_shift)
                    y1 = ybottom + (row * dY_spacing)
                    x2 = x1 + dX_spacing
                    y2 = y1
                    x3 = x1 + dX_shift
                    y3 = y1 + dY_spacing
                else:
                    x1 = xleft + (column * dX_shift)
                    y1 = ybottom + (row +1)* dY_spacing
                    x2 = x1 + dX_shift
                    y2 = y1 - dY_shift
                    x3 = x1 + dX_spacing
                    y3 = y1

            else:
                if row % 2 == 0:
                    x1 = xleft + column *  dX_shift
                    y1 = ybottom + (row + 1)* dY_spacing
                    x2 = x1 + dX_shift
                    y2 = y1 - dY_shift
                    x3 = x1 + dX_spacing
                    y3 = y1
                else:
                    x1 = xleft + column *  dX_shift
                    y1 = ybottom + (row )* dY_spacing
                    x2 = x1 + dX_spacing
                    y2 = y1
                    x3 = x1 + dX_shift
                    y3 = y1 + dY_spacing


            x = list()
            x.append(x1)
            x.append(x2)
            x.append(x3)

            y = list()
            y.append(y1)
            y.append(y2)
            y.append(y3)

            x_new , y_new = reproject_coordinates_batch(x, y, pProjection_reference_in)
            x1=x_new[0]
            x2=x_new[1]
            x3=x_new[2]

            y1=y_new[0]
            y2=y_new[1]
            y3=y_new[2]

            ring = ogr.Geometry(ogr.wkbLinearRing)
            ring.AddPoint(x1, y1)
            ring.AddPoint(x2, y2)
            ring.AddPoint(x3, y3)
            ring.AddPoint(x1, y1)
            pPolygon = ogr.Geometry(ogr.wkbPolygon)
            pPolygon.AddGeometry(ring)
            aCoords_gcs = np.full((4,2), -9999.0, dtype=float)
            aCoords_gcs[0,0] = x1
            aCoords_gcs[0,1] = y1
            aCoords_gcs[1,0] = x2
            aCoords_gcs[1,1] = y2
            aCoords_gcs[2,0] = x3
            aCoords_gcs[2,1] = y3
            aCoords_gcs[3,0] = x1
            aCoords_gcs[3,1] = y1


            dLongitude_center = float(np.mean(aCoords_gcs[0:3,0]))
            dLatitude_center = float(np.mean(aCoords_gcs[0:3,1]))

            iFlag = False
            if pPolygon.Within(pBoundary):
                iFlag = True
            else:
                #then check intersection
                if pPolygon.Intersects(pBoundary):
                    iFlag = True
                else:
                    pass

            if ( iFlag == True ):
                pTriangular = convert_gcs_coordinates_to_cell(6, dLongitude_center, dLatitude_center, aCoords_gcs)
                pTriangular.lCellID = lCellID
                dArea = pTriangular.calculate_cell_area()
                pTriangular.dArea = dArea
                pTriangular.calculate_edge_length()
                aTriangular.append(pTriangular)
                pFeature.SetGeometry(pPolygon)
                pFeature.SetField("cellid", lCellID)
                pFeature.SetField("longitude", dLongitude_center )
                pFeature.SetField("latitude", dLatitude_center )
                pFeature.SetField("area", dArea )
                pLayer.CreateFeature(pFeature)
                lCellID = lCellID + 1
            pass

    pDataset = pLayer = pFeature  = None


    return aTriangular
