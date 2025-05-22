
import geojson
import numpy as np
from osgeo import ogr
from jigsawpy import jigsaw_msh_t, savemsh


def linegeo(line, nset, eset, nobj, last):
    """
    LINEGEO: read a geoJSON line into a jigsaw msh_t object.

    """
    # Authors: Darren Engwirda

    npts = len(line) - 0

    if (npts > 0):
#----------------------------------------- read LINE dataset
        temp = jigsaw_msh_t()
        temp.vert2 = np.zeros(
            (npts + 0), dtype=temp.VERT2_t)
        temp.edge2 = np.zeros(
            (npts - 1), dtype=temp.EDGE2_t)

        temp.edge2["IDtag"][:] = nobj

        nobj = nobj + 1

        indx = np.arange(0, npts - 1) + last

        last = last + npts

        temp.vert2["coord"] = line[+0::]

        temp.edge2["index"][:, 0] = indx + 0
        temp.edge2["index"][:, 1] = indx + 1

        nset.append(temp.vert2)
        eset.append(temp.edge2)

    return nobj, last


def polygeo(loop, nset, eset, nobj, last):
    """
    POLYGEO: read a geoJSON poly into a jigsaw msh_t object.

    """
    # Authors: Darren Engwirda

    npts = len(loop) - 1

    if (npts > +0):
#----------------------------------------- read LOOP dataset
        temp = jigsaw_msh_t()
        temp.vert2 = np.zeros(
            (npts + 0), dtype=temp.VERT2_t)
        temp.edge2 = np.zeros(
            (npts - 0), dtype=temp.EDGE2_t)

        temp.edge2["IDtag"][:] = nobj

        nobj = nobj + 1

        indx = np.arange(0, npts - 1) + last

        itop = last + npts - 1

        idx1 = np.full(
            npts, itop, dtype=np.int32)
        idx1[:-1] = indx + 0

        idx2 = np.full(
            npts, last, dtype=np.int32)
        idx2[:-1] = indx + 1

        last = last + npts

        temp.vert2["coord"] = loop[:-1:]

        temp.edge2["index"][:, 0] = idx1
        temp.edge2["index"][:, 1] = idx2

        nset.append(temp.vert2)
        eset.append(temp.edge2)

    return nobj, last





def loadgeo(name, mesh):
    """
    LOADGEO: read a geoJSON file into a jigsaw msh_t object using GDAL.

    """

    nobj = +0; last = +0; nset = []; eset = []
    # Open the GeoJSON file
    driver = ogr.GetDriverByName("GeoJSON")
    datasource = driver.Open(name, 0)  # 0 means read-only

    if datasource is None:
        raise FileNotFoundError(f"Could not open {name}")

    # Get the layer from the datasource
    layer = datasource.GetLayer()

    for feature in layer:
        geom = feature.GetGeometryRef()

        if geom.GetGeometryType() == ogr.wkbLineString:
            line = geom.GetPoints()
            nobj, last = linegeo(line, nset, eset, nobj, last)

        elif geom.GetGeometryType() == ogr.wkbMultiLineString:
            for i in range(geom.GetGeometryCount()):
                line = geom.GetGeometryRef(i).GetPoints()
                nobj, last = linegeo(line, nset, eset, nobj, last)

        elif geom.GetGeometryType() == ogr.wkbPolygon:
            for i in range(geom.GetGeometryCount()):
                loop = geom.GetGeometryRef(i).GetPoints()
                nobj, last = polygeo(loop, nset, eset, nobj, last)

        elif geom.GetGeometryType() == ogr.wkbMultiPolygon:
            for i in range(geom.GetGeometryCount()):
                poly = geom.GetGeometryRef(i)
                for j in range(poly.GetGeometryCount()):
                    loop = poly.GetGeometryRef(j).GetPoints()
                    nobj, last = polygeo(loop, nset, eset, nobj, last)



    mesh.vert2 = np.concatenate(nset, axis=0)
    mesh.edge2 = np.concatenate(eset, axis=0)

    return

if (__name__ == "__main__"):


    mesh = jigsaw_msh_t()
    sFilename = "/qfs/people/liao313/data/hexwatershed/susquehanna/vector/hydrology/boundary_wgs.geojson"

    loadgeo(name=sFilename, mesh=mesh)


    sFilename = "test_new.msh"

    savemsh(name=sFilename, mesh=mesh)


