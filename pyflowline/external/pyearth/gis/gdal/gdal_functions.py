import os,sys
from math import cos, sin, sqrt, acos
import numpy as np
import osgeo
from osgeo import ogr, osr, gdal
#from shapely.wkt import loads
#most of these functions are copied from the pyearth package

import importlib
iFlag_cython = importlib.util.find_spec("cython") 
if iFlag_cython is not None:
    from pyflowline.algorithms.cython.kernel import longlat_to_3d
else:
    from pyflowline.algorithms.auxiliary.longlat_to_3d import longlat_to_3d



#https://stackoverflow.com/questions/8204998/how-to-check-if-a-pointlonc-latc-lie-on-a-great-circle-running-from-lona-lata
def calculate_distance_to_plane(x1, y1, x2, y2, x3, y3):
    
    a = np.radians(np.array((x1, y1) ))
    b = np.radians(np.array((x2, y2) )) #this is the middle one
    c = np.radians(np.array((x3, y3) ))
    # The points in 3D space
    a3 = longlat_to_3d(*a)
    b3 = longlat_to_3d(*b)
    c3 = longlat_to_3d(*c)
    #The formula is x+b*y+c*z=0 
    x1,y1,z1 = a3
    x2,y2,z2 = b3
    x3,y3,z3 = c3
    c = (-x1*y3 + x3* y1)/( z1*y3 - z3*y1 )
    b = (  -x1*z3 + x3 * z1 ) / (y1 * z3 - y3*z1)
    distance = abs(  x2 + b * y2 + c * z2 )
    return distance

def calculate_angle_betwen_vertex_normal(x1, y1, x2, y2, x3, y3):
    a = np.radians(np.array((x1, y1) ))
    b = np.radians(np.array((x2, y2) ))
    c = np.radians(np.array((x3, y3) ))
    # The points in 3D space
    a3 = longlat_to_3d(*a)
    b3 = longlat_to_3d(*b)
    c3 = longlat_to_3d(*c)
    a3vec = a3 - b3
    c3vec = c3 - b3 
    dot=np.dot(a3vec, c3vec)
    g = np.cross(a3vec, c3vec)
    det = np.dot(b3, g)
    angle = np.arctan2(det, dot)
    f = np.degrees(angle) 
    if f < 0:
        f = 360 + f
    
    return f

def calculate_angle_betwen_vertex(x1, y1, x2, y2, x3, y3):
    #all in degree
    a = np.radians(np.array((x1, y1) ))
    b = np.radians(np.array((x2, y2) ))
    c = np.radians(np.array((x3, y3) ))
    # The points in 3D space
    a3 = longlat_to_3d(*a)
    b3 = longlat_to_3d(*b)
    c3 = longlat_to_3d(*c)
    # Vectors in 3D space
    a3vec = a3 - b3
    c3vec = c3 - b3        
    angle3deg = angle_between_vectors_degrees(a3vec, c3vec)
    return  angle3deg
    
def angle_between_vectors_degrees(u, v):
    """Return the angle between two vectors in any dimension space,
    in degrees.
    """
    a = np.dot(u, v)
    b = np.linalg.norm(u)
    c = np.linalg.norm(v)
    d = a / (b* c)
    if d > 1:
        d = 1
    if d < -1:
        d = -1
    e = acos(d)
    f = np.degrees(e)
    return f
    
def convert_360_to_180(dLongitude_in):
    """[This function is modified from
    http://www.idlcoyote.com/map_tips/lonconvert.html]

    Args:
        dLongitude_in ([type]): [description]

    Returns:
        [type]: [description]
    """ 
    a = int(dLongitude_in /180)
    dLongitude_out = dLongitude_in - a*360.0
    return dLongitude_out

def convert_180_to_360(dLongitude_in):
    """[This function is modified from
    http://www.idlcoyote.com/map_tips/lonconvert.html]

    Args:
        dLongitude_in ([type]): [description]

    Returns:
        [type]: [description]
    """

    dLongitude_out = (dLongitude_in + 360.0) % 360.0

    return dLongitude_out

def retrieve_geotiff_metadata(sFilename_geotiff_in):
    """[retrieve the metadata of a geotiff file]

    Args:
        sFilename_geotiff_in ([type]): [description]

    Returns:
        [type]: [description]
    """
    pDriver = gdal.GetDriverByName('GTiff')
   
    pDataset = gdal.Open(sFilename_geotiff_in, gdal.GA_ReadOnly)

    if pDataset is None:
        print("Couldn't open this file: " + sFilename_geotiff_in)
        sys.exit("Try again!")
    else: 
        pProjection = pDataset.GetProjection()
        pSpatial_reference = osr.SpatialReference(wkt=pProjection)    
        ncolumn = pDataset.RasterXSize
        nrow = pDataset.RasterYSize        
        pGeotransform = pDataset.GetGeoTransform()
        dOriginX = pGeotransform[0]
        dOriginY = pGeotransform[3]
        dPixelWidth = pGeotransform[1]
        pPixelHeight = pGeotransform[5]               
        return dPixelWidth, dOriginX, dOriginY, nrow, ncolumn, pSpatial_reference, pProjection, pGeotransform

def degree_to_meter(dLatitude_in, dResolution_degree_in):
    """[Conver a degree-based resolution to meter-based resolution]

    Args:
        dLatitude_in ([type]): [description]
        dResolution_degree_in ([type]): [description]

    Returns:
        [type]: [description]
    """
    dRadius = 6378137.0  #unit: m earth radius
    dRadius2 = dRadius * np.cos( dLatitude_in / 180.0 * np.pi)
    dResolution_meter = dResolution_degree_in / 360.0 * (2*np.pi * dRadius2)

    return dResolution_meter

def meter_to_degree(dResolution_meter_in, dLatitude_mean_in):
    """[Convert a meter-based resolution to degree-based resolution]

    Args:
        dResolution_meter_in ([type]): [description]
        dLatitude_mean_in ([type]): [description]

    Returns:
        [type]: [description]
    """
    dLatitude_mean = abs(dLatitude_mean_in)

    dRadius = 6378137.0 # #unit: m earth radius
    dRadius2 = dRadius * np.cos( dLatitude_mean / 180.0 * np.pi)

    ##dResolution_meter = dResolution_degree / 360.0 * 2*np.pi * dRadius2

    dResolution_degree= dResolution_meter_in/(2*np.pi * dRadius2) * 360.0

    return dResolution_degree

def get_utm_spatial_reference(dLongitude_in):
    if -180 <= dLongitude_in <= 180:
        zone = int((dLongitude_in + 180) / 6) + 1
        hemisphere = 'N' if dLongitude_in >= 0 else 'S'
        epsg_code = 32600 + zone if hemisphere == 'N' else 32700 + zone
        utm_sr = osr.SpatialReference()
        utm_sr.ImportFromEPSG(epsg_code)
        return utm_sr
    else:
        raise ValueError("Longitude must be in the range [-180, 180].")

def reproject_coordinates(x_in, y_in, spatial_reference_source, spatial_reference_target=None):
    """[Reproject coordinates from one reference to another. By default to WGS84.]

    Args:
        x_in ([type]): [description]
        y_in ([type]): [description]
        spatial_reference_source ([type]): [description]
        spatial_reference_target ([type], optional): [description]. Defaults to None.

    Returns:
        [type]: [description]
    """      
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
   
    x_new,y_new, z = pTransform.TransformPoint( x_in,y_in)
    
    return x_new,y_new

def reproject_coordinates_batch(aX_in, aY_in, spatial_reference_source, spatial_reference_target=None):
    """[Reproject coordinates from one reference to another in batch mode. By default to WGS84.]

    Args:
        aX_in ([type]): [description]
        aY_in ([type]): [description]
        spatial_reference_source ([type]): [description]
        spatial_reference_target ([type], optional): [description]. Defaults to None.

    Returns:
        [type]: [description]
    """
    #Reproject a list of x,y coordinates. 

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

    npoint = len(aX_in)
    x_new=list()
    y_new=list()
    for i in range(npoint):
        x0 = aX_in[i]
        y0 = aY_in[i]
   
        x1,y1, z = pTransform.TransformPoint( x0,y0)

        x_new.append(x1)
        y_new.append(y1)
    
    return x_new,y_new

def calculate_distance_based_on_lon_lat(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * np.arcsin(sqrt(a)) 
    # Radius of earth in kilometers. Use 3956 for miles
    r = 6378137.0
    return c * r

def calculate_polygon_area(aLon_in, aLat_in,  algorithm = 0):
    """
    Computes area of spherical polygon, assuming spherical Earth. 
    Returns result in ratio of the sphere's area if the radius is specified. Otherwise, in the units of provided radius.
    lats and lons are in degrees.
    """
    #TODO: take into account geodesy (i.e. convert latitude to authalic sphere, use radius of authalic sphere instead of mean radius of spherical earth)
    #for i in range(len(aLon_in)):
    #    aLon_in[i] = aLon_in[i] + 180.0
    
    lons = np.deg2rad(aLon_in)
    lats = np.deg2rad(aLat_in)
    if algorithm==0:
        # Line integral based on Green's Theorem, assumes spherical Earth      
        #close polygon
        if lats[0]==lats[-1] and lons[0]==lons[-1] :
            pass
        else:
            lats = np.append(lats, lats[0])
            lons = np.append(lons, lons[0])

        # Get colatitude (a measure of surface distance as an angle)
        a = np.sin(lats/2)**2 + np.cos(lats)* np.sin(lons/2)**2
        colat = 2*np.arctan2( np.sqrt(a), np.sqrt(1-a) )

        #azimuth of each point in segment from the arbitrary origin
        az = np.arctan2(np.cos(lats) * np.sin(lons), np.sin(lats)) % (2*np.pi)

        # Calculate step sizes
        # daz = np.diff(az) % (2*pi)
        daz = np.diff(az)
        daz = (daz + np.pi) % (2 * np.pi) - np.pi

        # Determine average surface distance for each step
        deltas=np.diff(colat)/2
        colat=colat[0:-1]+deltas

        # Integral over azimuth is 1-cos(colatitudes)
        integrands = (1-np.cos(colat)) * daz

        # Integrate and save the answer as a fraction of the unit sphere.
        # Note that the sum of the integrands will include a factor of 4pi.
        area = abs(sum(integrands))/(4*np.pi) # Could be area of inside or outside

        area = min(area,1-area)
        radius= 6378137.0
       
        dArea  = area * 4*np.pi*(radius**2)
        return dArea


    elif algorithm==2:
        #L'Huilier Theorem, assumes spherical earth
        #see:
        # https://mathworld.wolfram.com/SphericalPolygon.html
        # https://web.archive.org/web/20160324191929/http://forum.worldwindcentral.com/showthread.php?20724-A-method-to-compute-the-area-of-a-spherical-polygon
        # https://github.com/spacetelescope/spherical_geometry/blob/master/spherical_geometry/polygon.py
        # https://github.com/tylerjereddy/spherical-SA-docker-demo/blob/master/docker_build/demonstration.py
        #TODO
        pass
    elif algorithm==3:
        #https://trs.jpl.nasa.gov/handle/2014/41271
        #TODO
        pass

def gdal_read_geotiff_file(sFilename_in):
    """Read a Geotiff format raster file.

    Args:
        sFilename_in (string): The file name

    Returns:
        tuple: aData_out, pPixelWidth, dOriginX, dOriginY, nrow, ncolumn, dMissing_value , pGeotransform, pProjection,  pSpatial_reference
    """
    
    if os.path.exists(sFilename_in):
        pass
    else:
        print('The file does not exist!')
        return

    sDriverName='GTiff'
    pDriver = gdal.GetDriverByName(sDriverName)  

    if pDriver is None:
        print ("%s pDriver not available.\n" % sDriverName)
    else:
        print  ("%s pDriver IS available.\n" % sDriverName)  

    pDataset = gdal.Open(sFilename_in, gdal.GA_ReadOnly)

    if pDataset is None:
        print("Couldn't open this file: " + sFilename_in)
        sys.exit("Try again!")
    else:       
        pProjection = pDataset.GetProjection()

        pDataset.GetMetadata()
       
        ncolumn = pDataset.RasterXSize
        nrow = pDataset.RasterYSize
        nband = pDataset.RasterCount

        pGeotransform = pDataset.GetGeoTransform()
        dOriginX = pGeotransform[0]
        dOriginY = pGeotransform[3]
        dPixelWidth = pGeotransform[1]
        pPixelHeight = pGeotransform[5]

        pBand = pDataset.GetRasterBand(1)

        # Data type of the values
        gdal.GetDataTypeName(pBand.DataType)
        # Compute statistics if needed
        if pBand.GetMinimum() is None or pBand.GetMaximum() is None:
            pBand.ComputeStatistics(0)

        dMissing_value = pBand.GetNoDataValue()
       
        aData_out = pBand.ReadAsArray(0, 0, ncolumn, nrow)
    
        #we will use one of them to keep the consistency
        pSpatial_reference = osr.SpatialReference(wkt=pProjection)
       

        pDataset = None
        pBand = None      
        pBand = None

        return aData_out, dPixelWidth, dOriginX, dOriginY, nrow, ncolumn, dMissing_value, pGeotransform, pProjection,  pSpatial_reference

def Google_MetersPerPixel( zoomLevel ):  
   
   # Return to the caller if there is an error.
   #On_Error, 2
   
   #; Need a zoom level?
   
   
   # Number of pixels in an image with a zoom level of 0.
   pixels_in_image = 256
   
   # The equitorial radius of the Earth assuming WGS-84 ellipsoid.
   earth_radius = 6378137.0
   
   # The number of meters per pixel.
   metersPerPixel = (2* np.pi * earth_radius) / pixels_in_image / np.power(2,zoomLevel)
   
   # Return the value.
   return metersPerPixel

def read_mesh_boundary(sFilename_boundary_in):
    """
    convert a shpefile to json format.
    This function should be used for stream flowline only.
    """
    iReturn_code = 1
    if os.path.isfile(sFilename_boundary_in):
        pass
    else:
        print('This mesh file does not exist: ', sFilename_boundary_in )
        iReturn_code = 0
        return iReturn_code

    
    pDriver_json = ogr.GetDriverByName('GeoJSON')    
    pDataset_mesh = pDriver_json.Open(sFilename_boundary_in, gdal.GA_ReadOnly)
    pLayer_mesh = pDataset_mesh.GetLayer(0)
    pSpatial_reference_out = pLayer_mesh.GetSpatialRef()
    ldefn = pLayer_mesh.GetLayerDefn()   

    #we also need to spatial reference
    for pFeature_mesh in pLayer_mesh:
        pGeometry_mesh = pFeature_mesh.GetGeometryRef()                     
        pGeometrytype_boundary = pGeometry_mesh.GetGeometryName()
        if(pGeometrytype_boundary == 'POLYGON'):       
            pBoundary_ogr = pGeometry_mesh  
        else:
            if(pGeometrytype_boundary == 'MULTIPOLYGON'):    
                nLine = pGeometry_mesh.GetGeometryCount()
                for i in range(nLine):
                    pBoundary_ogr = pGeometry_mesh.GetGeometryRef(i)
               
                pass
            else:
                pass
            pass   
            
            
    pBoundary_wkt = pBoundary_ogr.ExportToWkt()
    aExtent = pBoundary_ogr.GetEnvelope()
    min_x, max_x, min_y, max_y = aExtent
   
    return pBoundary_wkt, aExtent
    

def get_geometry_coords(geometry):
    
    sGeometry_type = geometry.GetGeometryName()
    if sGeometry_type =='POINT':
        return get_point_coords(geometry)
    elif sGeometry_type == 'LINESTRING':
        return get_linestring_coords(geometry)
    elif sGeometry_type =='POLYGON':
        return get_polygon_exterior_coords(geometry)
    else:
        raise ValueError("Unsupported geometry type.")

def get_polygon_exterior_coords(polygon_geometry):
    exterior_coords = []
    ring = polygon_geometry.GetGeometryRef(0)  # Get the exterior ring
    for i in range(ring.GetPointCount()):
        point = ring.GetPoint(i)
        exterior_coords.append((point[0], point[1]))
    return np.array(exterior_coords)

def get_linestring_coords(linestring_geometry):
    coords = []
    for i in range(linestring_geometry.GetPointCount()):
        point = linestring_geometry.GetPoint(i)
        coords.append((point[0], point[1]))
    return np.array(coords)

def get_point_coords(point_geometry):
    point = point_geometry.GetPoint()
    return np.array([(point[0], point[1])])