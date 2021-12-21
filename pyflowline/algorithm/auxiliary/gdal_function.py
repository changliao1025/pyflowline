import sys
from osgeo import gdal, osr
import numpy as np
def obtain_raster_metadata(sFilename_geotiff):
    pDriver = gdal.GetDriverByName('GTiff')
   
    pDataset = gdal.Open(sFilename_geotiff, gdal.GA_ReadOnly)

    if pDataset is None:
        print("Couldn't open this file: " + sFilename_geotiff)
        sys.exit("Try again!")
    else: 
        pProjection = pDataset.GetProjection()
        pSpatialRef = osr.SpatialReference(wkt=pProjection)
    
    
        ncolumn = pDataset.RasterXSize
        nrow = pDataset.RasterYSize
        #nband = pDataset.RasterCount

        pGeotransform = pDataset.GetGeoTransform()
        dOriginX = pGeotransform[0]
        dOriginY = pGeotransform[3]
        dPixelWidth = pGeotransform[1]
        pPixelHeight = pGeotransform[5]       
        
        print( dPixelWidth, dOriginX, dOriginY, nrow, ncolumn)
        return dPixelWidth, dOriginX, dOriginY, nrow, ncolumn, pSpatialRef, pProjection, pGeotransform


def degree_to_meter(dLatitude, dResolution_degree):
    dRadius = 6378.1 * 1000 #unit: m earth radius
    dRadius2 = dRadius * np.cos( dLatitude / 180.0 * np.pi)
    dResolution_meter = dResolution_degree / 360.0 * (2*np.pi * dRadius2)

    return dResolution_meter

def meter_to_degree(dResolution_meter, dLatitude_mean):
    dLatitude_mean = abs(dLatitude_mean)

    dRadius = 6378.1 * 1000 #unit: m earth radius
    dRadius2 = dRadius * np.cos( dLatitude_mean / 180.0 * np.pi)

    ##dResolution_meter = dResolution_degree / 360.0 * 2*np.pi * dRadius2

    dResolution_degree= dResolution_meter/(2*np.pi * dRadius2) * 360.0

    return dResolution_degree