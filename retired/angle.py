import sys
from telnetlib import X3PAD
import numpy as np
import osgeo
import math
from math import radians, cos, sin, asin, sqrt
from osgeo import ogr, osr, gdal, gdalconst
from numpy import arctan2, cos, sin, sqrt, pi, power, append, diff

def calculate_angle_betwen_vertex_2d_new(x1, y1, x2, y2, x3, y3):
    a = np.radians(np.array((x1, y1) ))
    b = np.radians(np.array((x2, y2) ))
    c = np.radians(np.array((x3, y3) ))
    # The points in 3D space
    a3 = longlat_to_3d(*a)
    b3 = longlat_to_3d(*b)
    c3 = longlat_to_3d(*c)

    a3vec = a3 - b3
    c3vec = c3 - b3 

    x1,y1,z1=a3vec
    x2,y2,z2=c3vec



    xn, yn, zn = b3
    dot = x1*x2 + y1*y2 + z1*z2
    h0 = np.dot(a3vec, c3vec)
    print(dot-h0)

    det = x1*y2*zn + x2*yn*z1 + xn*y1*z2 - z1*y2*xn - z2*yn*x1 - zn*y1*x2

    g = np.cross(a3vec, c3vec)
    det1 = np.dot(b3, g)
    print(det-det1)

    angle = np.arctan2(det, dot)
    f = np.degrees(angle) 
    
    
    return f

def calculate_angle_betwen_vertex_2d(x1, y1, x2, y2, x3, y3):
    a0 = np.arctan2(y3 - y2, x3 - x2)
    a00 = np.degrees(a0)
    b0 = np.arctan2(y1 - y2, x1 - x2)
    b00 = np.degrees(b0)
    angle = ( a0 - b0  )
    angle0 = a00-b00
    
    

    a = np.radians(np.array((x1, y1) ))
    b = np.radians(np.array((x2, y2) ))
    c = np.radians(np.array((x3, y3) ))

    # Vectors in latitude/longitude space
    avec = a - b
    cvec = c - b

    # Adjust vectors for changed longitude scale at given latitude into 2D space
    lat = b[0]
    avec[1] *= math.cos(lat)
    cvec[1] *= math.cos(lat)
    # Find the angle between the vectors in 2D space
    angle2deg = angle_between_vectors_degrees(avec, cvec)

    if angle < 0:
        c = np.degrees(angle)
        d = c +  360
        e = 360 - angle2deg
    else:
        pass
    
    return angle2deg
def angle_between(c1, c2, c3):
    c1 = np.radians(np.array(c1 ))
    c2 = np.radians(np.array(c2 ))
    c3 = np.radians(np.array(c3 ))
    
    angle = (np.arctan2(c3[1] - c1[1], c3[0] - c1[0]) -
             np.arctan2(c2[1] - c1[1], c2[0] - c1[0]))
    
    if angle < 0:
        angle +=  2 * np.pi
    
    angle = np.degrees(angle)
    return angle

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
    
def longlat_to_3d(lonr, latr):
    """Convert a point given latitude and longitude in radians to
    3-dimensional space, assuming a sphere radius of one."""
    a = math.cos(latr) * math.cos(lonr)
    b = math.cos(latr) * math.sin(lonr)
    c = math.sin(latr)
    return np.array((a,b,c))

def longlat_to_3d2(lonr, latr):
    """Convert a point given latitude and longitude in radians to
    3-dimensional space, assuming a sphere radius of one."""
    a = math.cos(latr) * math.cos(lonr)
    b = math.cos(latr) * math.sin(lonr)
    c = math.sin(latr)
    return np.array((b,a,c))

def angle_between_vectors_degrees(u, v):
    """Return the angle between two vectors in any dimension space,
    in degrees."""

    a = np.dot(u, v)
    b = np.linalg.norm(u)
    c = np.linalg.norm(v)
    d = a / (b * c)
    if d > 1:
        d = 1
    if d < -1:
        d = -1
    e = math.acos(d)
    f = np.degrees(e) 

    return f

x1=-76.6633525247768
y1=40.0755204044615
x2=-76.65805228803946
y2=40.071942697076814
x3=-76.6544421247906
y3=40.0695058044708

a = calculate_angle_betwen_vertex_2d_new(x1, y1, x2, y2, x3, y3)
print(a)

c= calculate_angle_betwen_vertex(x1, y1, x2, y2, x3, y3)
print(c)

#b= angle_between([x2,y2],[x1,y1],[x3,y3])
#print(b)