import os
import numpy as np
from osgeo import  osr, gdal, ogr
#from shapely.wkt import loads
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import matplotlib.ticker as mticker
import matplotlib.path as mpath
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from pyflowline.external.pyearth.toolbox.math.stat.remap import remap
from pyflowline.external.pyearth.gis.gdal.gdal_functions import get_geometry_coords
class OOMFormatter(mpl.ticker.ScalarFormatter):
    def __init__(self, order=0, fformat="%1.1e", offset=True, mathText=True):
        self.oom = order
        self.fformat = fformat
        mpl.ticker.ScalarFormatter.__init__(self,useOffset=offset,useMathText=mathText)
    def _set_order_of_magnitude(self):
        self.orderOfMagnitude = self.oom
    def _set_format(self, vmin=None, vmax=None):
        self.format = self.fformat
        if self._useMathText:
            self.format = r'$ mathdefault{%s}$' % self.format

def map_vector_polyline_data(sFilename_in,
                             sFilename_output_in= None,
                             iFlag_color_in = None,
                             iFlag_label_in = None,
                             iFlag_thickness_in =None,
                             iFont_size_in = None,
                             sField_thickness_in=None,
                             sField_color_in=None,
                             iFlag_scientific_notation_colorbar_in=None,
                             sColormap_in = None,
                             sTitle_in = None,
                             iDPI_in = None,
                             dMissing_value_in=None,
                             dData_max_in = None,
                             dData_min_in = None,
                             sFont_in = None,
                             sExtend_in =None,
                             sUnit_in=None,
                             aLegend_in = None,
                             aExtent_in = None,
                             pProjection_map_in = None):
    """
    Map a vector polyline data

    Args:
        iFiletype_in (int): _description_
        sFilename_in (_type_): _description_
        sFilename_output_in (_type_): _description_
        iFlag_thickness_in (_type_, optional): _description_. Defaults to None.
        sField_thickness_in (_type_, optional): _description_. Defaults to None.
        iFlag_scientific_notation_colorbar_in (_type_, optional): _description_. Defaults to None.
        sColormap_in (_type_, optional): _description_. Defaults to None.
        sTitle_in (_type_, optional): _description_. Defaults to None.
        iDPI_in (_type_, optional): _description_. Defaults to None.
        dMissing_value_in (_type_, optional): _description_. Defaults to None.
        dData_max_in (_type_, optional): _description_. Defaults to None.
        dData_min_in (_type_, optional): _description_. Defaults to None.
        sFont_in (_type_, optional): _description_. Defaults to None.
        sExtend_in (_type_, optional): _description_. Defaults to None.
        sUnit_in (_type_, optional): _description_. Defaults to None.
        aLegend_in (_type_, optional): _description_. Defaults to None.
        aExtent_in (_type_, optional): _description_. Defaults to None.
        pProjection_map_in (_type_, optional): _description_. Defaults to None.
    """

    if os.path.isfile(sFilename_in):
        pass
    else:
        print('The file does not exist: ', sFilename_in)
        return

    sFilename_out= sFilename_output_in
    if iDPI_in is not None:
        iDPI = iDPI_in
    else:
        iDPI = 300

    if iFlag_scientific_notation_colorbar_in is not None:
        iFlag_scientific_notation_colorbar = iFlag_scientific_notation_colorbar_in
    else:
        iFlag_scientific_notation_colorbar = 0

    if iFlag_color_in is not None:
        iFlag_color = iFlag_color_in
    else:
        iFlag_color = 0

    if iFlag_thickness_in is not None:
        iFlag_thickness = iFlag_thickness_in
    else:
        iFlag_thickness = 0

    if iFont_size_in is not None:
        iFont_size = iFont_size_in
    else:
        iFont_size = 12

    if dData_min_in is not None:
        iFlag_data_min = 1
        dData_min = dData_min_in
    else:
        iFlag_data_min = 0
        pass

    if dData_max_in is not None:
        iFlag_data_max = 1
        dData_max = dData_max_in
    else:
        iFlag_data_max = 0
        pass

    if sField_thickness_in is not None:
        sField_thickness = sField_thickness_in
    else:
        sField_thickness = ''

    pDriver = ogr.GetDriverByName('GeoJSON')

    pDataset = pDriver.Open(sFilename_in, gdal.GA_ReadOnly)
    pLayer = pDataset.GetLayer(0)


    if sColormap_in is not None:
        sColormap = sColormap_in
    else:
        sColormap =  'Spectral'

    if sTitle_in is not None:
        sTitle = sTitle_in
        iFlag_title =1
    else:
        iFlag_title=0
        sTitle =  ''

    if sExtend_in is not None:
        sExtend = sExtend_in
    else:
        sExtend =  'max'

    if sUnit_in is not None:
        sUnit = sUnit_in
    else:
        sUnit =  ''

    if sFont_in is not None:
        sFont = sFont_in
    else:
        sFont = "Times New Roman"

    plt.rcParams["font.family"] = sFont

    cmap = cm.get_cmap(sColormap)

    pSrs = osr.SpatialReference()
    pSrs.ImportFromEPSG(4326)    # WGS84 lat/lon

    lID = 0
    dLat_min = 90
    dLat_max = -90
    dLon_min = 180
    dLon_max = -180
    for pFeature in pLayer:
        pGeometry_in = pFeature.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()

        if sGeometry_type =='LINESTRING':
            #dummy0 = loads( pGeometry_in.ExportToWkt() )
            #aCoords_gcs = dummy0.coords
            #aCoords_gcs= np.array(aCoords_gcs)
            aCoords_gcs = get_geometry_coords(pGeometry_in)
            aCoords_gcs = aCoords_gcs[:,0:2]

            dLon_max = np.max( [dLon_max, np.max(aCoords_gcs[:,0])] )
            dLon_min = np.min( [dLon_min, np.min(aCoords_gcs[:,0])] )
            dLat_max = np.max( [dLat_max, np.max(aCoords_gcs[:,1])] )
            dLat_min = np.min( [dLat_min, np.min(aCoords_gcs[:,1])] )

    if pProjection_map_in is not None:
        pProjection_map = pProjection_map_in
    else:
        pProjection_map = ccrs.Orthographic(central_longitude =  0.50*(dLon_max+dLon_min),
                                            central_latitude = 0.50*(dLat_max+dLat_min),
                                            globe=None)


    fig = plt.figure( dpi = iDPI  )
    iSize_x= 8
    iSize_y= 8
    fig.set_figwidth( iSize_x )
    fig.set_figheight( iSize_y )
    ax = fig.add_axes([0.1, 0.15, 0.75, 0.8] , projection=pProjection_map ) #request.crs
    ax.set_global()

    nColor = pLayer.GetFeatureCount()
    aColor = cm.rainbow(np.linspace(0, 1, nColor))

    if iFlag_thickness ==1:
        aValue =list()
        for pFeature in pLayer:
            dValue = pFeature.GetField(sField_thickness)
            aValue.append(dValue)
            pass

        aValue = np.array(aValue)
        dValue_max = np.max(aValue)
        dValue_min = np.min(aValue)
        iThickness_max = 2.0
        iThickness_min = 0.3

    for pFeature in pLayer:
        pGeometry_in = pFeature.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()
        if iFlag_thickness ==1:
            dValue = pFeature.GetField(sField_thickness)

        if sGeometry_type =='LINESTRING':
            #dummy0 = loads( pGeometry_in.ExportToWkt() )
            #aCoords_gcs = dummy0.coords
            #aCoords_gcs= np.array(aCoords_gcs)
            aCoords_gcs = get_geometry_coords(pGeometry_in)
            aCoords_gcs = aCoords_gcs[:,0:2]
            nvertex = len(aCoords_gcs)
            if nvertex == 2 :
                dLon_label = 0.5 * (aCoords_gcs[0][0] + aCoords_gcs[1][0] )
                dLat_label = 0.5 * (aCoords_gcs[0][1] + aCoords_gcs[1][1] )
            else:
                lIndex_mid = int(nvertex/2)
                dLon_label = aCoords_gcs[lIndex_mid][0]
                dLat_label = aCoords_gcs[lIndex_mid][1]

            codes = np.full(nvertex, mpath.Path.LINETO, dtype=int )
            codes[0] = mpath.Path.MOVETO
            path = mpath.Path(aCoords_gcs, codes)
            x, y = zip(*path.vertices)

            if iFlag_thickness ==1:
                iThickness = remap( dValue, dValue_min, dValue_max, iThickness_min, iThickness_max )
            else:
                iThickness = 1.0

            if iFlag_color ==1:
                if nColor < 10:
                    rgba = aColor[lID]
                else:
                    iColor_index = (dValue-dValue_min ) /(dValue_max - dValue_min )
                    rgba = cmap(iColor_index)

            else:
                rgba='black'

            line, = ax.plot(x, y, color=rgba, linewidth=iThickness, transform=ccrs.Geodetic())
            lID = lID + 1

    if aExtent_in is None:
        marginx  = (dLon_max - dLon_min) / 50
        marginy  = (dLat_max - dLat_min) / 50
        aExtent = [dLon_min - marginx , dLon_max + marginx , dLat_min - marginy , dLat_max + marginy]
    else:
        aExtent = aExtent_in

    ax.set_extent(aExtent)
    ax.coastlines(color='black', linewidth=1)
  
    if aLegend_in is not None:
        nlegend = len(aLegend_in)
        dLocation0 = 0.96
        for i in range(nlegend):
            sText = aLegend_in[i]
            dLocation = dLocation0 - i * 0.06
            ax.text(0.03, dLocation, sText,
                    verticalalignment='top', horizontalalignment='left',
                    transform=ax.transAxes,
                    color='black',
                    fontsize=iFont_size )

            pass

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=0.5, linestyle='--')

    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlocator = mticker.MaxNLocator(5)
    gl.ylocator = mticker.MaxNLocator(5)
    gl.xlabel_style = {'size': 8, 'color': 'k', 'rotation':0, 'ha':'right'}
    gl.ylabel_style = {'size': 8, 'color': 'k', 'rotation':90,'weight': 'normal'}


    if iFlag_title==1:
        ax.set_title( sTitle )
    else:
        pass


    pDataset = pLayer = pFeature  = None
    if sFilename_output_in is None:
        plt.show()
    else:
        sDirname = os.path.dirname(sFilename_output_in)
        sFilename = os.path.basename(sFilename_output_in)
        sFilename_out = os.path.join(sDirname, sFilename)
        sExtension = os.path.splitext(sFilename)[1]
        if sExtension == '.png':
            plt.savefig(sFilename_out, bbox_inches='tight')
        else:
            if sExtension == '.pdf':
                plt.savefig(sFilename_out, bbox_inches='tight')
            else:
                plt.savefig(sFilename_out, bbox_inches='tight', format ='ps')
    #clean
    plt.close('all')
    plt.clf()
