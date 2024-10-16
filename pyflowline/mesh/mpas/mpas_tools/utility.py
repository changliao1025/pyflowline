
import numpy as np
import jigsawpy

from pyflowline.mesh.mpas.jigsaw.inpoly2 import inpoly2


def addpoly(geom, poly, itag):
    """
    ADDPOLY: add new closed polygon POLY to mst_t obj. GEOM.

    The POLY.POINT + POLY.EDGE2 arrays are appended to GEOM,
    and a new "loop" added to GEOM.BOUND. This new loop is
    assigned ID = ITAG.

    """
    # Authors: Darren Engwirda

    temp = jigsawpy.jigsaw_msh_t()

    temp.point = poly.point
    temp.edge2 = poly.edge2

    zipmesh(temp)                       # ensure compressed

    temp.bound = np.empty(
        (poly.edge2.size), dtype=geom.BOUND_t)
    temp.bound["index"] = \
        np.arange(0, poly.edge2.size)
    temp.bound["cells"] = \
        jigsawpy.jigsaw_def_t.JIGSAW_EDGE2_TAG
    temp.bound["IDtag"] = itag

    temp.edge2["index"] += geom.vert2.size
    temp.bound["index"] += geom.edge2.size

    geom.point = np.concatenate(
        (geom.vert2, temp.vert2), axis=0)
    geom.edge2 = np.concatenate(
        (geom.edge2, temp.edge2), axis=0)
    geom.bound = np.concatenate(
        (geom.bound, temp.bound), axis=0)

    return


def addline(geom, line, itag):
    """
    ADDLINE: push new open polyline LINE to mst_t obj. GEOM.

    The LINE.POINT + LINE.EDGE2 arrays are appended to GEOM.
    The new edges are assigned ID = ITAG.

    """
    # Authors: Darren Engwirda

    temp = jigsawpy.jigsaw_msh_t()

    temp.point = line.point
    temp.edge2 = line.edge2

    temp.edge2["IDtag"] = itag

    zipmesh(temp)                       # ensure compressed

    temp.edge2["index"] += geom.vert2.size

    geom.point = np.concatenate(
        (geom.vert2, temp.vert2), axis=0)
    geom.edge2 = np.concatenate(
        (geom.edge2, temp.edge2), axis=0)

    return


def zipmesh(mesh):
    """
    ZIPMESH: "zip" a mst_t obj., pruning any unused points /
    cells, and compressing cell indexing.

    """
    # Authors: Darren Engwirda

    used = np.full(
        mesh.point.size, False, dtype=bool)

#---------------------------------- flag nodes used in cells
    if (mesh.edge2 is not None and
            mesh.edge2.size > +0):

        used[mesh.edge2[
            "index"].reshape(-1)] = True

    if (mesh.tria3 is not None and
            mesh.tria3.size > +0):

        used[mesh.tria3[
            "index"].reshape(-1)] = True

    if (mesh.quad4 is not None and
            mesh.quad4.size > +0):

        used[mesh.quad4[
            "index"].reshape(-1)] = True

#---------------------------------- re-index cells, compress
    redo = np.full(
        mesh.point.size, 0, dtype=np.int32)

    redo[used] = np.arange(
        0, np.count_nonzero(used))

    if (mesh.edge2 is not None and
            mesh.edge2.size > +0):

        mesh.edge2["index"] = \
            redo[mesh.edge2["index"]]

    if (mesh.tria3 is not None and
            mesh.tria3.size > +0):

        mesh.tria3["index"] = \
            redo[mesh.tria3["index"]]

    if (mesh.quad4 is not None and
            mesh.quad4.size > +0):

        mesh.quad4["index"] = \
            redo[mesh.quad4["index"]]

#---------------------------------- prune any un-used points
    mesh.point = mesh.point[used]

    return


def innerto(vert, geom):
    """
    INNERTO: return ID-tags of polygons enclosing each point
    in VERT. Return -1 for exterior points.

    """
    # Authors: Darren Engwirda

    itag = np.full(
        vert.shape[+0], -1, dtype=np.int32)

    imin = np.min(geom.bound["IDtag"])
    imax = np.max(geom.bound["IDtag"])

    for ipos in range(imin, imax + 1):

#---------------------------------- inpolygon for POLY=ITAG
        indx = geom.bound["index"][
            geom.bound["IDtag"] == ipos]

        xpts = geom.vert2["coord"]
        edge = geom.edge2["index"][indx]

        mask, _ = inpoly2(vert, xpts, edge)

        itag[mask] = ipos

    return itag
