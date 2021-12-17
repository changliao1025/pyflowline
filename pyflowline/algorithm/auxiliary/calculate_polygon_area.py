import numpy as np
from numpy import arctan2, cos, sin, sqrt, pi, power, append, diff
def calculate_polygon_area(lats, lons, algorithm = 0, radius = 6378137.0):
    """
    Computes area of spherical polygon, assuming spherical Earth. 
    Returns result in ratio of the sphere's area if the radius is specified. Otherwise, in the units of provided radius.
    lats and lons are in degrees.
    """
    #TODO: take into account geodesy (i.e. convert latitude to authalic sphere, use radius of authalic sphere instead of mean radius of spherical earth)
    lats = np.deg2rad(lats)
    lons = np.deg2rad(lons)

    if algorithm==0:
        # Line integral based on Green's Theorem, assumes spherical Earth
        

        #close polygon
        if lats[0]!=lats[-1]:
            lats = append(lats, lats[0])
            lons = append(lons, lons[0])

        # Get colatitude (a measure of surface distance as an angle)
        a = sin(lats/2)**2 + cos(lats)* sin(lons/2)**2
        colat = 2*arctan2( sqrt(a), sqrt(1-a) )

        #azimuth of each point in segment from the arbitrary origin
        az = arctan2(cos(lats) * sin(lons), sin(lats)) % (2*pi)

        # Calculate step sizes
        # daz = diff(az) % (2*pi)
        daz = diff(az)
        daz = (daz + pi) % (2 * pi) - pi

        # Determine average surface distance for each step
        deltas=diff(colat)/2
        colat=colat[0:-1]+deltas

        # Integral over azimuth is 1-cos(colatitudes)
        integrands = (1-cos(colat)) * daz

        # Integrate and save the answer as a fraction of the unit sphere.
        # Note that the sum of the integrands will include a factor of 4pi.
        area = abs(sum(integrands))/(4*pi) # Could be area of inside or outside

        area = min(area,1-area)
        if radius is not None: #return in units of radius
            return area * 4*pi*radius**2
        else: #return in ratio of sphere total area
            return area
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