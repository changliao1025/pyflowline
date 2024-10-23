
import numpy as np
from scipy import spatial
from scipy import interpolate
import jigsawpy

from pyflowline.mesh.mpas.jigsaw.utility import zipmesh


def sphdist(rsph, xone, yone, xtwo, ytwo):
    """
    SPHDIST: return the distance from the points [XONE,YONE]
    to the point [XONE,YONE] on a sphere of radii RSPH.

    """
    # Authors: Darren Engwirda

    dlon = .5 * (xone - xtwo)
    dlat = .5 * (yone - ytwo)

    dist = 2. * rsph * np.arcsin(np.sqrt(
        np.sin(dlat) ** 2 +
        np.sin(dlon) ** 2 * np.cos(yone) * np.cos(ytwo)
    ))

    return dist


def rd_dist(xone, xtwo):
    """
    RD_DIST: return the distance between the points XONE and
    XTWO in R^D. Points are N-by-D arrays.

    """
    # Authors: Darren Engwirda

    return np.sqrt(np.sum((xtwo - xone) ** 2, axis=1,
                          keepdims=True))


def blender(val1, val2, dist, blen, bgap):
    """
    BLENDER: return a blend of VAL1 and VAL2, as a nonlinear
    weighted combination.

    TANH((DIST-LEN) ./ GAP) weighting is used -- blending is
    centred about DIST=LEN and spans an approx. width of GAP.

    """
    # Authors: Darren Engwirda

    beta = .5 + .5 * np.tanh(
        np.pi * (dist - blen) / bgap)

    return (+0.0 + beta) * val2 + (+1.0 - beta) * val1


def R3toS2(radii, E3):
    """
    R3TOS2: return the LON-LAT coord's associated with a set
    of points in R^3. A (geocentric) projection is first
    done to ensure the points lie on the ellipsoidal surface,
    with the projected points then transformed to [LON,LAT]
    pairs.

    """
    # Authors: Darren Engwirda

    PP = .5 * E3

    ax = PP[:, 0] ** 1 / radii[0] ** 1
    ay = PP[:, 1] ** 1 / radii[1] ** 1
    az = PP[:, 2] ** 1 / radii[2] ** 1

    aa = ax ** 2 + ay ** 2 + az ** 2

    bx = PP[:, 0] ** 2 / radii[0] ** 2
    by = PP[:, 1] ** 2 / radii[1] ** 2
    bz = PP[:, 2] ** 2 / radii[2] ** 2

    bb = bx * 2. + by * 2. + bz * 2.

    cx = PP[:, 0] ** 1 / radii[0] ** 1
    cy = PP[:, 1] ** 1 / radii[1] ** 1
    cz = PP[:, 2] ** 1 / radii[2] ** 1

    cc = cx ** 2 + cy ** 2 + cz ** 2
    cc = cc - 1.0

    ts = bb * bb - 4. * aa * cc

    ok = ts >= .0

    AA = aa[ok]; BB = bb[ok]; CC = cc[ok]; TS = ts[ok]

    t1 = (-BB + np.sqrt(TS)) / AA / 2.0
    t2 = (-BB - np.sqrt(TS)) / AA / 2.0

    tt = np.maximum(t1, t2)

    P3 = np.zeros(E3.shape, dtype=float)
    P3[ok, 0] = (1. + tt) * PP[ok, 0]
    P3[ok, 1] = (1. + tt) * PP[ok, 1]
    P3[ok, 2] = (1. + tt) * PP[ok, 2]

    return jigsawpy.R3toS2(radii, P3)


def zipnear(init, geom, spac, near=10.0):
    """
    ZIPNEAR: "zip" nodes / cells in INIT. that are "too near"
    to geometry features. Here, "too near" means they lie
    within NEAR * H(x) of objects in GEOM. where H(x) is the
    spacing distribution defined via SPAC.

    The INIT. object is modified in-place.

    """
    # Authors: Darren Engwirda

#------------------------------------ subdiv. geom. for H(x)

    divgeom(geom, spac)

#------------------------------------ find dist. to geometry

    tree = spatial.cKDTree(
        jigsawpy.S2toR3(
            geom.radii, geom.point["coord"]))

    dmax = np.max(spac.value) * near

    dist, _ = tree.query(
        init.point["coord"],
        eps=0.0, distance_upper_bound=dmax)

    apos = jigsawpy.R3toS2(
        geom.radii, init.point["coord"][:])

#------------------------------------ zip init. if too close

    hfun = interpolate.RectBivariateSpline(
        spac.ygrid, spac.xgrid, spac.value)

    hval = hfun(
        apos[:, 1], apos[:, 0], grid=False)

    keep = dist >= hval * near

#------------------------------------ re-index init. for zip

    if (init.edge2 is not None and
            init.edge2.size > +0):

        mask = np.logical_and.reduce((
            keep[init.edge2["index"][:, 0]],
            keep[init.edge2["index"][:, 1]]
        ))

        init.edge2 = init.edge2[mask]

    if (init.tria3 is not None and
            init.tria3.size > +0):

        mask = np.logical_and.reduce((
            keep[init.tria3["index"][:, 0]],
            keep[init.tria3["index"][:, 1]],
            keep[init.tria3["index"][:, 2]]
        ))

        init.tria3 = init.tria3[mask]

    zipmesh(init)

    return


def divgeom(geom, spac):
    """
    DIVGEOM: subdivide the edges in the msh_t object GEOM to
    satisfy the mesh spacing definition given by SPAC.

    The GEOM. object is modified in-place.

    """
    # Authors: Darren Engwirda

    kind = geom.mshID.lower()

    hfun = interpolate.RectBivariateSpline(
        spac.ygrid, spac.xgrid, spac.value)

    while True:                             # while too long

        vert = geom.point["coord"]
        cell = geom.edge2["index"]

    #-------------------------------- eval. spacing at nodes

        hval = hfun(
            vert[:, 1], vert[:, 0], grid=False)

    #-------------------------------- eval. edge length, x^d

        if (kind == "ellipsoid-mesh"):

            xpos = jigsawpy.S2toR3(
                geom.radii, geom.point["coord"])

            elen = sphdist(
                np.mean(geom.radii),
                vert[cell[:, 0], 0],
                vert[cell[:, 0], 1],
                vert[cell[:, 1], 0],
                vert[cell[:, 1], 1])

        if (kind == "euclidean-mesh"):

            xpos = np.array(vert[:], copy=True)

            elen = rd_dist(
                vert[cell[:, 0], :],
                vert[cell[:, 1], :])

        hmid = 0.5 * (
            hval[cell[:, 0]] + hval[cell[:, 1]]
        )

        mask = elen >= hmid * 4. / 3.       # TRUE to subdiv

        ndiv = np.count_nonzero(mask)

        if (not np.any(mask)): break

        xmid = 0.5 * (
            xpos[cell[:, 0]] + xpos[cell[:, 1]]
        )

    #-------------------------------- subdiv. edge at middle

        if (kind == "ellipsoid-mesh"):

            newx = np.zeros(
                (ndiv), dtype=geom.point.dtype)

            newx["coord"] = R3toS2(
                geom.radii, xmid[mask])

        if (kind == "euclidean-mesh"):

            newx = np.zeros(
                (ndiv), dtype=geom.point.dtype)

            newx["coord"] = xmid[mask]

    #-------------------------------- re-index for new edges

        inew = np.arange(
            +0, ndiv) + geom.point.size

        new1 = np.empty(
            (ndiv), dtype=geom.EDGE2_t)
        new2 = np.empty(
            (ndiv), dtype=geom.EDGE2_t)

        new1["index"][:, 0] = \
            geom.edge2["index"][mask, 0]
        new1["index"][:, 1] = inew
        new1["IDtag"] = \
            geom.edge2["IDtag"][mask]

        new2["index"][:, 0] = inew
        new2["index"][:, 1] = \
            geom.edge2["index"][mask, 1]
        new2["IDtag"] = \
            geom.edge2["IDtag"][mask]

        geom.edge2[mask] = new1

    #-------------------------------- add new cells to GEOM.

        geom.point = np.concatenate(
            (geom.point, newx), axis=0)
        geom.edge2 = np.concatenate(
            (geom.edge2, new2), axis=0)

    return
