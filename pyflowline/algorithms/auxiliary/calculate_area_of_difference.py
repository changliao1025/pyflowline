import os
import json
import numpy as np
from osgeo import ogr, osr, gdal, gdalconst
from shapely.geometry import Point, LineString, MultiLineString
from shapely.wkt import loads

from shapely.ops import polygonize, polygonize_full

from pyflowline.classes.edge import pyedge
from pyflowline.classes.vertex import pyvertex
from pyflowline.classes.flowline import pyflowline

from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_flowline 

from pyflowline.algorithms.intersect.intersect_flowline_with_flowline import intersect_flowline_with_flowline
from pyflowline.algorithms.auxiliary.gdal_functions import calculate_polygon_area

def calculate_area_of_difference_raw(sFilename_a, sFilename_b):

    #not yet supported


    return

def calculate_area_of_difference_simplified(aFlowline_in, \
     sFilename_output_in):

    if os.path.exists(sFilename_output_in): 
        #delete it if it exists
        os.remove(sFilename_output_in)
        pass

    pDriver_geojson = ogr.GetDriverByName( "GeoJSON")

    nFlowline = len(aFlowline_in)
    aFlowline = list()
    
    for i in range(nFlowline):
       
        pFlowline = aFlowline_in[i]
        nVertex= pFlowline.nVertex
        aCoords_gcs = np.full((nVertex,2), -9999. ,dtype=float)
        for k in range(nVertex):
            aCoords_gcs[k,0] = pFlowline.aVertex[k].dLongitude_degree
            aCoords_gcs[k,1] = pFlowline.aVertex[k].dLatitude_degree
        
        aCoords_gcs= tuple(j for j in aCoords_gcs) #np.array(aCoords_gcs)    
        aFlowline.append(aCoords_gcs)
    
    

    dummy = polygonize(aFlowline)
    aPolygon_out = list(dummy)
    pDataset = pDriver_geojson.CreateDataSource(sFilename_output_in) 
    pSpatial_reference_gcs = osr.SpatialReference()  
    pSpatial_reference_gcs.ImportFromEPSG(4326)    # WGS84 lat/lon
    pSpatial_reference_gcs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

    pLayer = pDataset.CreateLayer('intersect', pSpatial_reference_gcs, ogr.wkbPolygon)
    # Add one attribute
    pLayer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger64)) #long type for high resolution

    pLayerDefn = pLayer.GetLayerDefn()
    pFeature = ogr.Feature(pLayerDefn)    
    lCellID =0
    dArea =0.0

    
    for po in aPolygon_out:
        ring = ogr.Geometry(ogr.wkbLinearRing)
        aCoords_gcs = po.exterior.coords
        aCoords_gcs= np.array(aCoords_gcs)  
        nPoint  = (aCoords_gcs.shape)[0]
        lons=list()
        lats=list()
        for i in range(nPoint):
            ring.AddPoint(aCoords_gcs[i,0], aCoords_gcs[i,1])
            lons.append( aCoords_gcs[i,0] )
            lats.append( aCoords_gcs[i,1] )
        
        pPolygon = ogr.Geometry(ogr.wkbPolygon)
        pPolygon.AddGeometry(ring)
        pFeature.SetGeometry(pPolygon)
        pFeature.SetField("id", lCellID)
        pLayer.CreateFeature(pFeature)
        lCellID= lCellID+1

        dArea0 = calculate_polygon_area(lons, lats)
        dArea = dArea + dArea0

    return aPolygon_out, dArea
   
    

    

      

             

    
