import numpy as np
from scipy.spatial import ConvexHull

def find_minimal_enclosing_polygon(aLongitude_degree, aLatitude_degree):

    # Convert aLatitude_degree/aLongitude_degree points to Cartesian coordinates
    aVertex_cartesian = list()
    for i in range(len(aLongitude_degree)):
        #x = np.cos(np.deg2rad(aLatitude_degree[i])) * np.cos(np.deg2rad(aLongitude_degree[i]))
        #y = np.cos(np.deg2rad(aLatitude_degree[i])) * np.sin(np.deg2rad(aLongitude_degree[i]))
        #z = np.sin(np.deg2rad(aLatitude_degree[i]))
        #aVertex_cartesian.append((x,y,z))
        x = aLongitude_degree[i]
        y = aLatitude_degree[i]
        aVertex_cartesian.append((x,y))


    aVertex_cartesian = np.array(aVertex_cartesian)
    # Create a ConvexHull object from the Cartesian coordinates
    hull = ConvexHull(aVertex_cartesian)
    # Get the indices of the points on the convex hull
    indices = hull.vertices
    # Get the points on the convex hull in Cartesian coordinates
    aVertex_on_hull_cartesian = [aVertex_cartesian[i] for i in indices]
    # Convert the points on the convex hull back to aLatitude_degree/aLongitude_degree coordinates
    aVertex_on_hull_latlon = []
    for p in aVertex_on_hull_cartesian:
        #dLongitude_degree = np.rad2deg(np.arctan2(p[1], p[0]))
        #dLatitude_degree = np.rad2deg(np.arcsin(p[2]))
        dLongitude_degree = p[0]
        dLatitude_degree = p[1]
        aVertex_on_hull_latlon.append((dLongitude_degree, dLatitude_degree))

    return aVertex_on_hull_latlon

