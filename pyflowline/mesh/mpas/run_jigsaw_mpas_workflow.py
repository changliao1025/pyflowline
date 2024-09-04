import os
from pathlib import Path
import math
import copy
import time
import numpy as np
import jigsawpy
from pyflowline.mesh.mpas.inpoly2 import inpoly2
from pyflowline.mesh.mpas.saveesm import saveesm
def run_jigsaw_mpas_workflow(sWorkspace_jigsaw_in,
                             aConfig_in = None,
            projector=[0.0, 0.0],
            iFlag_geom=True,
              iFlag_spac=True,
                iFlag_init=True,
                iFlag_opts=True,
                  aFilename_river_in=None,
                               aFilename_watershed_boundary_in= None,
                               aFilenamae_lake_boundary_in = None,
                               aFilename_coastline_in = None):
    """
    run_jigsaw_mpas_workflow: main call to the ICoM mesh-gen. infrastructure.

    1. Call user-defined MESH-PATH/COMPOSE.py to assemble
       mesh geometry, spacing pattern, initial conditions
       and mesh-gen. parameters.
    2. Call JIGSAW to build the triangulation.
    4. Call MPAS meshtools to make the MPAS mesh data files.

    """
    # Authors: Darren Engwirda and Chang Liao

    class obj(object): pass                 # dummy object

    aFlag_option = obj()
    aFlag_option.geom = iFlag_geom
    aFlag_option.init = iFlag_init
    aFlag_option.spac = iFlag_spac
    aFlag_option.opts = iFlag_opts

#---------- call JIGSAW to build the initial triangular mesh

    geom, gprj, mesh, mprj =  runjgsw(sWorkspace_jigsaw_in, aFlag_option, projector, aConfig_in=aConfig_in)

#---------- call utilities to write output for ATS and FVCOM

    #if (isinstance(idtag_ats, list) and len(idtag_ats) > 0       and idtag_ats[0] >= +0):

#-------------------------------------- write output for ATS
        #mout = zipmesh(mprj, idtag_ats)
        #saveats(sWorkspace_jigsaw_in, mout)

    #if (isinstance(idtag_fvc, list) and len(idtag_fvc) > 0       and idtag_fvc[0] >= +0):

#-------------------------------------- write output for FVC
        #mout = zipmesh(mprj, idtag_fvc)
    #   savefvc(sWorkspace_jigsaw_in, mout)

#---------- run MPAS meshtools to build MPAS data-structures



#-------------------------------------- write output for ESM
    saveesm(sWorkspace_jigsaw_in, geom, mesh)



    #jigsawpy.cmd.jigsaw(opts, mesh)

    return

def runjgsw(sWorkspace_jigsaw_in, aFlag_option, projector, aConfig_in=None):
    """
    RUNJGSW: main call to JIGSAW to build the triangulation.

    MESH-PATH should point to a user-defined mesh directory,
    containing the COMPOSE.py template.

    Firstly, MESH-PATH/COMPOSE.py is called to build
    user-defined geometry, initial conditions, mesh spacing
    and mesh configuration information.

    The boolean flags MAKE-BOOL control whether mesh
    information is built from scratch, or if an exitsing an
    existing file is to be used. For example, setting
    MAKE-BOOL.SPAC = FALSE relies on an existing spacing
    pattern be available in MESH-PATH/tmp/.

    This information is written to MESH-PATH/tmp/ to be
    accessed by subsequent calls to JIGSAW.

    Finally, JIGSAW is run to build the triangular mesh,
    calling either the multi-level (TETRIS) or single-level
    (JIGSAW) algorithms. Cells are assigned ID-tags via
    the polygon/regions defined in GEOM.BOUNDS.

    Returns full-dimensional and 2d-projected msh_t objects.

    """
    # Authors: Darren Engwirda

    mesh = jigsawpy.jigsaw_msh_t()
    mprj = jigsawpy.jigsaw_msh_t()
    gprj = jigsawpy.jigsaw_msh_t()

#------------------------------------ setup via user COMPOSE
    FULL_SPHERE_RADIUS = +6371.0 #unit check

    if (aFlag_option.geom):
        #geom = getattr(import_module(
        #    base + ".compose"), "setgeom")()
        geom = jigsawpy.jigsaw_msh_t()
        if aConfig_in is not None:
            geom.mshID = aConfig_in["geom_mshID"]
            geom.radii = np.full(3, aConfig_in['FULL_SPHERE_RADIUS'], dtype=geom.REALS_t)
        else:
            geom.mshID = "ellipsoid-mesh"
            geom.radii = np.full(3, FULL_SPHERE_RADIUS, dtype=geom.REALS_t)


    if (aFlag_option.spac):
        #spac = getattr(import_module(
        #    base + ".compose"), "setspac")()
        spac = jigsawpy.jigsaw_msh_t()
        if aConfig_in is not None:
            spac.mshID = aConfig_in["spac_mshID"]
            spac.radii = np.full(  3, aConfig_in['FULL_SPHERE_RADIUS'], dtype=spac.REALS_t)
            spac.xgrid = np.linspace( -1. * np.pi, +1. * np.pi, 360)
            spac.ygrid = np.linspace( -.5 * np.pi, +.5 * np.pi, 180)
            spac.value = np.full(  (180, 360), +5.0E+001, dtype=spac.REALS_t)
        else:
            spac.mshID = "ellipsoid-grid"
            spac.radii = np.full( 3, aConfig_in['FULL_SPHERE_RADIUS'], dtype=spac.REALS_t)
            spac.xgrid = np.linspace(  -1. * np.pi, +1. * np.pi, 360)
            spac.ygrid = np.linspace( -.5 * np.pi, +.5 * np.pi, 180)
            spac.value = np.full(  (180, 360), +5.0E+001, dtype=spac.REALS_t)


    if (aFlag_option.init):
        #init = getattr(import_module(
        #    base + ".compose"), "setinit")()
        init = jigsawpy.jigsaw_msh_t()

    if (aFlag_option.opts):
        #opts = getattr(import_module(
        #    base + ".compose"), "setopts")()
        opts = jigsawpy.jigsaw_jig_t()
        if aConfig_in is not None:
            opts.hfun_scal = aConfig_in["hfun_scal"]
            opts.hfun_hmax = aConfig_in["hfun_hmax"]
            opts.hfun_hmin = aConfig_in["hfun_hmin"]
            opts.mesh_dims = aConfig_in["mesh_dims"]
            opts.bisection = aConfig_in["bisection"]
            opts.optm_qlim = aConfig_in["optm_qlim"]
            opts.optm_iter = aConfig_in["optm_iter"]
            opts.optm_qtol = aConfig_in["optm_qtol"]
        else:
            opts.hfun_scal = "absolute"
            opts.hfun_hmax = float("inf")       # null spacing lim
            opts.hfun_hmin = float(+0.00)
            opts.mesh_dims = +2                 # 2-dim. simplexes
            opts.bisection = -1                 # call heutristic!
            opts.optm_qlim = +9.5E-01           # tighter opt. tol
            opts.optm_iter = +32
            opts.optm_qtol = +1.0E-05

#------------------------------------ setup files for JIGSAW

    sWorkspace_tmp = os.path.join(
        sWorkspace_jigsaw_in, "tmp")
    Path(sWorkspace_tmp).mkdir(parents=True, exist_ok=True)

    opts.geom_file = os.path.join(
        sWorkspace_jigsaw_in, "tmp", "geom.msh")

    opts.jcfg_file = os.path.join(
        sWorkspace_jigsaw_in, "tmp", "opts.jig")

    opts.init_file = os.path.join(
        sWorkspace_jigsaw_in, "tmp", "init.msh")

    opts.hfun_file = os.path.join(
        sWorkspace_jigsaw_in, "tmp", "spac.msh")

    opts.mesh_file = os.path.join(
        sWorkspace_jigsaw_in, "tmp", "mesh.msh")

    opts.hfun_tags = "precision = 9"    # less float prec.

    jigsawpy.savemsh(opts.geom_file, geom,
                     opts.geom_tags)

    jigsawpy.savemsh(opts.hfun_file, spac,
                     opts.hfun_tags)

    jigsawpy.savemsh(opts.init_file, init,
                     opts.init_tags)

#------------------------------------ make mesh using JIGSAW

    if (not hasattr(opts, "bisection")):

        opts.bisection = +0

    if (opts.bisection < +0):           # bisect heuristic

        rbar = np.mean(geom.radii)
        hbar = np.mean(spac.value)

        nlev = round(math.log2(
            rbar / math.sin(.4 * math.pi) / hbar)
        )

        nlev = nlev - 1

        ttic = time.time()

        jigsawpy.cmd.tetris(opts, nlev - 0, mesh)

        ttoc = time.time()

        print("CPUSEC =", (ttoc - ttic))
        print("BISECT =", +nlev)

    elif (opts.bisection > +0):         # bisect specified

        nlev = opts.bisection

        ttic = time.time()

        jigsawpy.cmd.tetris(opts, nlev - 0, mesh)

        ttoc = time.time()

        print("CPUSEC =", (ttoc - ttic))
        print("BISECT =", +nlev)

    else:                               # do non-recursive

        ttic = time.time()

        jigsawpy.cmd.jigsaw(opts, mesh)

        ttoc = time.time()

        print("CPUSEC =", (ttoc - ttic))

#------------------------------------ form local projections

    gprj = copy.deepcopy(geom)          # local 2d objects
    mprj = copy.deepcopy(mesh)

    if (mesh.vert3.size > +0):
        project(geom, mesh, gprj, mprj, projector)

#------------------------------------ assign IDtag's to cell

    if (geom.bound is not None and
            geom.bound.size > +0):      # tags per polygon

        imin = np.amin(geom.bound["IDtag"])
        imax = np.amax(geom.bound["IDtag"])

        for itag in range(
                imin + 0, imax + 1):

            tagcell(geom, mesh, gprj, mprj, itag)

#------------------------------------ check mesh for quality

    cost = jigsawpy.triscr2(            # quality metrics!
        mesh.point["coord"],
        mesh.tria3["index"])

    print("TRISCR =", np.min(cost), np.mean(cost))

    cost = jigsawpy.pwrscr2(
        mesh.point["coord"],
        mesh.power,
        mesh.tria3["index"])

    print("PWRSCR =", np.min(cost), np.mean(cost))

    tbad = jigsawpy.centre2(
        mesh.point["coord"],
        mesh.power,
        mesh.tria3["index"])

    print("OBTUSE =",
          +np.count_nonzero(np.logical_not(tbad)))

    ndeg = jigsawpy.trideg2(
        mesh.point["coord"],
        mesh.tria3["index"])

    print("TOPOL. =",
          +np.count_nonzero(ndeg==+6) / ndeg.size)

#------------------------------------ save mesh for Paraview

    sWorkspace_tmp = os.path.join(
        sWorkspace_jigsaw_in, "out")
    Path(sWorkspace_tmp).mkdir(parents=True, exist_ok=True)

    jigsawpy.savevtk(os.path.join(
        sWorkspace_jigsaw_in, "out", "geom.vtk"), geom)
    jigsawpy.savevtk(os.path.join(
        sWorkspace_jigsaw_in, "out", "spac.vtk"), spac)
    jigsawpy.savevtk(os.path.join(
        sWorkspace_jigsaw_in, "out", "init.vtk"), init)
    jigsawpy.savevtk(os.path.join(
        sWorkspace_jigsaw_in, "out", "mesh.vtk"), mesh)

    jigsawpy.savevtk(os.path.join(
        sWorkspace_jigsaw_in, "out", "geom_prj.vtk"), gprj)
    jigsawpy.savevtk(os.path.join(
        sWorkspace_jigsaw_in, "out", "mesh_prj.vtk"), mprj)

    return geom, gprj, mesh, mprj


def project(geom, mesh, gprj, mprj, pmid):
    """
    PROJECT: projection of GEOM/MESH objects to a 2d. plane.

    Modifies GPRJ, MPRJ objects "inplace".

    """
    # Authors: Darren Engwirda

    proj = jigsawpy.jigsaw_prj_t()

    mprj.point = np.full(   mesh.point.size, +0, dtype=mprj.VERT2_t)

    mprj.point["coord"] =   jigsawpy.R3toS2(
            geom.radii, mesh.point["coord"])

    proj.prjID = "stereographic"
    proj.radii = np.mean(geom.radii)
    proj.xbase = pmid[0] * np.pi / +180.
    proj.ybase = pmid[1] * np.pi / +180.

    jigsawpy.project(gprj, proj, "fwd")
    jigsawpy.project(mprj, proj, "fwd")

    return


def tagcell(geom, mesh, gprj, mprj, itag):
    """
    TAGCELL: assign ID-tags to mesh cells based-on polygons
    defined in GEOM.BOUND.

    Modifies MESH, MPRJ objects "inplace".

    """
    # Authors: Darren Engwirda

    this = geom.bound["IDtag"] == itag
    cell = geom.bound["index"][this]

    loop = geom.edge2["index"][cell, :]

#------------------------------------ compute cell "centres"
    ball = jigsawpy.tribal2(
        mprj.point["coord"], mesh.tria3["index"])

#------------------------------------ calc. in-polygon tests
    tmsk, _ = inpoly2(
        ball[:, :-1], gprj.point["coord"], loop)

    vmsk, _ = inpoly2(mprj.point["coord"],
                      gprj.point["coord"], loop)

    vbnd = np.full(
        mesh.point.size, False, dtype=np.bool)

    vbnd[mesh.edge2["index"].flatten()] = True

    vmsk[vbnd] = False

    tbnd = np.logical_and.reduce((      # all vert. on bnd
        vbnd[mesh.tria3["index"][:, 0]],
        vbnd[mesh.tria3["index"][:, 1]],
        vbnd[mesh.tria3["index"][:, 2]]))

    tmsk[np.logical_not(tbnd)] = False

    mask = np.logical_or.reduce((       # vert. or cell in
        vmsk[mesh.tria3["index"][:, 0]],
        vmsk[mesh.tria3["index"][:, 1]],
        vmsk[mesh.tria3["index"][:, 2]],
        tmsk))

#------------------------------------ assign IDtags to cells
    mesh.tria3["IDtag"][mask == 1] = itag
    mprj.tria3["IDtag"][mask == 1] = itag

    return


def zipmesh(mesh, tags):
    """
    ZIPMESH: "zip" a mesh down to just the verts/edges/cells
    associated with the list of ID's in TAGS.

    Returns a new "zipped" msh_t object.

    """
    # Authors: Darren Engwirda

    mout = jigsawpy.jigsaw_msh_t()

    keep_cell = np.full(
        mesh.tria3.size, False, dtype=np.bool)
    keep_edge = np.full(
        mesh.edge2.size, False, dtype=np.bool)
    keep_vert = np.full(
        mesh.point.size, False, dtype=np.bool)

#------------------------------------ "flag" tagged entities
    for itag in tags:

        keep_cell[mesh.tria3["IDtag"] == itag] = True

    keep_vert[mesh.tria3[
        "index"][keep_cell].flatten()] = True

    keep_edge = np.logical_and.reduce((
        keep_vert[mesh.edge2["index"][:, 0]],
        keep_vert[mesh.edge2["index"][:, 1]])
    )

    mout.point = mesh.point[keep_vert]
    mout.edge2 = mesh.edge2[keep_edge]
    mout.tria3 = mesh.tria3[keep_cell]

#------------------------------------ update vertex indexing
    redo = \
        np.zeros(mesh.point.size, dtype=np.int32)

    redo[keep_vert] = \
        np.arange(0, np.count_nonzero(keep_vert))

    mout.edge2["index"] = redo[mout.edge2["index"]]
    mout.tria3["index"] = redo[mout.tria3["index"]]

    return mout