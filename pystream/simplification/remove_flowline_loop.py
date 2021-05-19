import sys, os
import numpy as np
from osgeo import gdal, osr, ogr, gdalconst

from hexwatershed.preprocess.stream.simplification.add_unique_line import  add_unique_line
from hexwatershed.preprocess.stream.check_same_point  import check_same_point

def remove_flowline_loop(aFlowline_in):
    aFlowline_out=list()

    nFlowline = len(aFlowline_in)
        
    j = 0
    aLine =[]
    def find_paralle_stream(i, pt_start):
        iFlag =0
        ndownstream=0
        aDownstream=[]
        for j in range(nFlowline):
            pFeature_in = pLayer_in.GetFeature(j)            
            pGeometry_in = pFeature_in.GetGeometryRef()     
            npt = pGeometry_in.GetPointCount()            
            pt1_start = pGeometry_in.GetPoint(0)         
            pt1_end = pGeometry_in.GetPoint(npt-1)  
            if check_same_point(pt_start, pt1_start)==1 and j !=i :
                ndownstream= ndownstream+1
                aDownstream.append(j)
                
        return ndownstream, aDownstream
    lID=0
    aFlag = np.full(nFlowline, 0, dtype=int)
    for i in range(nFlowline):
        pFeature_in = pLayer_in.GetFeature(i)            
        pGeometry_in = pFeature_in.GetGeometryRef()     
        npt = pGeometry_in.GetPointCount()        
        pt_start = pGeometry_in.GetPoint(0)         
        pt_end = pGeometry_in.GetPoint(npt-1)  
        ndownstream , aDownstream = find_paralle_stream(i, pt_start)
        if ndownstream == 0:
            if aFlag[i] !=1:
                pFeatureOut.SetGeometry(pGeometry_in)
                pFeatureOut.SetField("id", lID)
                pLayer_out.CreateFeature(pFeatureOut)
                lID = lID + 1
                aFlag[i]=1
            pass
        else:                     
            #more than one, so we only take one
            if(ndownstream>0):
                
                if aFlag[i] !=1:
                    pFeatureOut.SetGeometry(pGeometry_in)
                    pFeatureOut.SetField("id", lID)
                    pLayer_out.CreateFeature(pFeatureOut)
                    lID = lID + 1
                    aFlag[i]=1
                for k in range(ndownstream):
                    aFlag[  aDownstream[k]] =1
            pass
        
        pass 

    return aFlowline_out