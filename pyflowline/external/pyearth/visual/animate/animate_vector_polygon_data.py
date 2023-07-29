import os
import json
import numpy as np
import cartopy.crs as ccrs
import matplotlib as mpl
#from shapely.wkt import loads
from osgeo import osr, gdal, ogr
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
from matplotlib import animation
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from pyflowline.external.pyearth.gis.gdal.gdal_functions import get_geometry_coords

pProjection = ccrs.PlateCarree()  # for latlon data only


class OOMFormatter(mpl.ticker.ScalarFormatter):
    def __init__(self, order=0, fformat="%1.1e", offset=True, mathText=True):
        self.oom = order
        self.fformat = fformat
        mpl.ticker.ScalarFormatter.__init__(
            self, useOffset=offset, useMathText=mathText
        )

    def _set_order_of_magnitude(self):
        self.orderOfMagnitude = self.oom

    def _set_format(self, vmin=None, vmax=None):
        self.format = self.fformat
        if self._useMathText:
            self.format = r"$\mathdefault{%s}$" % self.format


iFigwidth_default = 12
iFigheight_default = 12


def animate_vector_polygon_data(
    sFilename_mesh_in,
    sFilename_animation_json_in,
    sFilename_animation_out,
    iFlag_type_in=None,
    iFigwidth_in=None,
    iFigheight_in=None,
    aExtent_in=None,
    pProjection_map_in=None,
):
    if iFigwidth_in is None:
        iFigwidth_in = iFigwidth_default

    if iFigheight_in is None:
        iFigheight_in = iFigheight_default

    # read domain json
    # read result json

    # read animation json

    dLat_min = 90
    dLat_max = -90
    dLon_min = 180
    dLon_max = -180
    aCellID = list()
    aData = list()
    aCell_raw = list()
    aCell_new = list()
    aCell_animation = list()
    pDriver = ogr.GetDriverByName("GeoJSON")

    pDataset = pDriver.Open(sFilename_mesh_in, gdal.GA_ReadOnly)
    pLayer = pDataset.GetLayer(0)
    pSrs = osr.SpatialReference()
    pSrs.ImportFromEPSG(4326)  # WGS84 lat/lon
    for pFeature in pLayer:
        pGeometry_in = pFeature.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()

        if sGeometry_type == "POLYGON":
            #dummy0 = loads(pGeometry_in.ExportToWkt())
            #aCoords_gcs = dummy0.exterior.coords
            #aCoords_gcs = np.array(aCoords_gcs)
            aCoords_gcs = get_geometry_coords(pGeometry_in)

            dLon_max = np.max([dLon_max, np.max(aCoords_gcs[:, 0])])
            dLon_min = np.min([dLon_min, np.min(aCoords_gcs[:, 0])])
            dLat_max = np.max([dLat_max, np.max(aCoords_gcs[:, 1])])
            dLat_min = np.min([dLat_min, np.min(aCoords_gcs[:, 1])])

    if pProjection_map_in is not None:
        pProjection_map = pProjection_map_in
    else:
        pProjection_map = ccrs.Orthographic(
            central_longitude=0.50 * (dLon_max + dLon_min),
            central_latitude=0.50 * (dLat_max + dLat_min),
            globe=None,
        )

    fig = plt.figure(dpi=600)
    fig.set_figwidth(iFigwidth_in)
    fig.set_figheight(iFigheight_in)
    ax = fig.add_axes([0.1, 0.1, 0.65, 0.8], projection=pProjection_map)
    ax.set_global()

    aData = np.array(aData)
    dData_max = np.max(aData)
    dData_min = np.min(aData)
    marginx = (dLon_max - dLon_min) / 20
    marginy = (dLat_max - dLat_min) / 20
    if aExtent_in is None:
        aExtent = [
            dLon_min - marginx,
            dLon_max + marginx,
            dLat_min - marginy,
            dLat_max + marginy,
        ]
    else:
        aExtent = aExtent_in

    ax.set_extent(aExtent)
    ax.coastlines()  # resolution='110m')
    if iFlag_type_in == 1:  # full
        gl = ax.gridlines(
            crs=ccrs.PlateCarree(),
            draw_labels=True,
            linewidth=1,
            color="gray",
            alpha=0.3,
            linestyle="--",
        )
        gl.xlabel_style = {"size": 8, "color": "k", "rotation": 0, "ha": "right"}
        gl.ylabel_style = {"size": 8, "color": "k", "rotation": 90, "weight": "normal"}
    else:  # track mode
        pass

    cmap = cm.get_cmap("Spectral")
    # setting a title for the plot
    sText = "Priority flood in HexWatershed"
    ax.text(
        0.5,
        1.05,
        sText,
        verticalalignment="center",
        horizontalalignment="center",
        transform=ax.transAxes,
        color="black",
        fontsize=9,
    )
    cmap_reversed = cmap.reversed()

    with open(sFilename_mesh_in) as json_file:
        aCell_new = json.load(json_file)
        ncell_new = len(aCell_new)

    with open(sFilename_animation_json_in) as json_file:
        aCell_animation = json.load(json_file)
        ncell_animation = len(aCell_animation)

    sText = ""
    x1 = 0.0
    y1 = 0.0
    pArtist1 = ax.text(
        x1,
        y1,
        sText,
        verticalalignment="center",
        horizontalalignment="right",
        transform=ax.transAxes,
        color="black",
        fontsize=12,
    )

    x2 = 0.0
    y2 = 0.0
    pArtist2 = ax.text(
        x2,
        y2,
        sText,
        verticalalignment="center",
        horizontalalignment="left",
        transform=ax.transAxes,
        color="black",
        fontsize=12,
    )

    # trasform elevation
    norm = plt.Normalize(dData_min, dData_max)
    sm = plt.cm.ScalarMappable(cmap = cmap_reversed, norm=norm)
    sm.set_array(aData)
    ax_cb = fig.add_axes([0.85, 0.12, 0.04, 0.75])
    cb = fig.colorbar(sm, cax=ax_cb)
    sUnit = r"Unit: m"
    cb.ax.get_yaxis().set_ticks_position("right")
    cb.ax.get_yaxis().labelpad = 5
    cb.ax.set_ylabel(sUnit, rotation=90)
    cb.ax.get_yaxis().set_label_position("left")
    cb.ax.tick_params(labelsize=6)

    def animate(i):
        aArtist = list()
        # get the time step global id and updated elevation
        pcell_animation = aCell_animation[i]
        lCellID = int(pcell_animation["lCellID"])
        lCellIndex = aCellID.index(lCellID)
        pcell_new = aCell_new[lCellIndex]
        dlon = float(pcell_new["dLongitude_center_degree"])
        dlat = float(pcell_new["dLatitude_center_degree"])
        dummy0 = float(pcell_new["Elevation_raw"])
        dummy = float(pcell_new["Elevation"])
        avertex = pcell_new["vVertex"]
        nvertex = len(avertex)
        aLocation = np.full((nvertex, 2), 0.0, dtype=float)
        for k in range(nvertex):
            aLocation[k, 0] = avertex[k]["dLongitude_degree"]
            aLocation[k, 1] = avertex[k]["dLatitude_degree"]

        color_index = (dummy - dData_min) / (dData_max - dData_min)
        rgb = cmap_reversed(color_index)
        polygon = mpatches.Polygon(
            aLocation,
            closed=True,
            facecolor=rgb,
            edgecolor="none",
            transform=ccrs.Geodetic(),
        )
        pArtist0 = ax.add_patch(polygon)

        if iFlag_type_in == 1:
            dLon_min_zoom = dLon_min
            dLon_max_zoom = dLon_max
            dLat_min_zoom = dLat_min
            dLat_max_zoom = dLat_max
        else:
            dLon_min_zoom = dlon - marginx * 2
            dLon_max_zoom = dlon + marginx * 2
            dLat_min_zoom = dlat - marginy * 2
            dLat_max_zoom = dlat + marginy * 2

        if dummy0 > dummy:
            x1 = (dlon - dLon_min_zoom) / (dLon_max_zoom - dLon_min_zoom)
            y1 = (dlat - dLat_min_zoom) / (dLat_max_zoom - dLat_min_zoom) + 0.05
            x2 = (dlon - dLon_min_zoom) / (dLon_max_zoom - dLon_min_zoom)
            y2 = (dlat - dLat_min_zoom) / (dLat_max_zoom - dLat_min_zoom) - 0.05
        else:
            x1 = (dlon - dLon_min_zoom) / (dLon_max_zoom - dLon_min_zoom)
            y1 = (dlat - dLat_min_zoom) / (dLat_max_zoom - dLat_min_zoom) - 0.05
            x2 = (dlon - dLon_min_zoom) / (dLon_max_zoom - dLon_min_zoom)
            y2 = (dlat - dLat_min_zoom) / (dLat_max_zoom - dLat_min_zoom) + 0.05

        sText1 = "Before: " + "{:0.2f}".format(dummy0) + "m"

        pArtist1.set_x(x1)
        pArtist1.set_y(y1)
        pArtist1.set_text(sText1)
        sText2 = "After: " + "{:0.2f}".format(dummy) + "m"

        pArtist2.set_x(x2)
        pArtist2.set_y(y2)
        pArtist2.set_text(sText2)

        if iFlag_type_in == 2:
            aExtent_zoom = [
                dlon - marginx * 2,
                dlon + marginx * 2,
                dlat - marginy * 2,
                dlat + marginy * 2,
            ]
            ax.set_extent(aExtent_zoom)

        return pArtist0, pArtist1, pArtist2

    plt.rcParams[
        "animation.convert_path"
    ] = "/share/apps/ImageMagick/7.1.0-52/bin/convert"

    anim = animation.FuncAnimation(fig, animate, frames=ncell_animation, interval=400, blit=False)
    anim.save(sFilename_animation_out, writer="imagemagick")

    return
