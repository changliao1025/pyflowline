from osgeo import gdal, ogr
def check_file_type(filename):
    # Try to open as raster
    raster = gdal.Open(filename)
    if raster is not None:
        return 'raster'

    # Try to open as vector
    vector = ogr.Open(filename)
    if vector is not None:
        return 'vector'

    # If neither worked, the file is either not a spatial data file, or it's not a format that GDAL/OGR can read
    return 'unknown'

