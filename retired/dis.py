import sys
import numpy as np
import osgeo
import math
from math import radians, cos, sin, asin, sqrt
from osgeo import ogr, osr, gdal, gdalconst
from numpy import arctan2, cos, sin, sqrt, pi, power, append, diff

def calculate_distance_based_on_lon_lat(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    # Radius of earth in kilometers. Use 3956 for miles
    r = 6378137.0
    return c * r


x1 = -77.3881972192102
y1= 40.8147180155353

x2 = -77.0167967242282
y2= 40.3961176039638

x3 = -76.9213012661597
y3= 40.2893996309322

a =calculate_distance_based_on_lon_lat(x1, y1, x2, y2)
b =calculate_distance_based_on_lon_lat(x2, y2, x3, y3)
c =calculate_distance_based_on_lon_lat(x1, y1, x3, y3)

print(  c-(a+b))