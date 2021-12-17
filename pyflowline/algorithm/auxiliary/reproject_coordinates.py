import os, sys
import numpy as np
import osgeo
from osgeo import ogr, osr, gdal, gdalconst

def reproject_coordinates(x, y, spatial_reference_source, spatial_reference_target=None):
    """ Reproject a list of x,y coordinates. """

    if spatial_reference_target is not None:

        pass
    else:
        spatial_reference_target = osr.SpatialReference()
        spatial_reference_target.ImportFromEPSG(4326)
        
        pass

    
    if int(osgeo.__version__[0]) >= 3:
    # GDAL 3 changes axis order: https://github.com/OSGeo/gdal/issues/1546
                    
        spatial_reference_source.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)
        spatial_reference_target.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)

    
    pTransform = osr.CoordinateTransformation( spatial_reference_source, spatial_reference_target)
   
    x_new,y_new, z = pTransform.TransformPoint( x,y)
    
    return x_new,y_new

def reproject_coordinates_batch(x, y, spatial_reference_source, spatial_reference_target=None):
    """ Reproject a list of x,y coordinates. """

    if spatial_reference_target is not None:

        pass
    else:
        spatial_reference_target = osr.SpatialReference()
        spatial_reference_target.ImportFromEPSG(4326)
        
        pass

    
    if int(osgeo.__version__[0]) >= 3:
    # GDAL 3 changes axis order: https://github.com/OSGeo/gdal/issues/1546
                    
        spatial_reference_source.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)
        spatial_reference_target.SetAxisMappingStrategy(osgeo.osr.OAMS_TRADITIONAL_GIS_ORDER)

    
    pTransform = osr.CoordinateTransformation( spatial_reference_source, spatial_reference_target)

    npoint = len(x)
    x_new=list()
    y_new=list()
    for i in range(npoint):
        x0 = x[i]
        y0 = y[i]
   
        x1,y1, z = pTransform.TransformPoint( x0,y0)

        x_new.append(x1)
        y_new.append(y1)
    
    return x_new,y_new