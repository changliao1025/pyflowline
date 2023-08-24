import os
import numpy as np
from osgeo import  osr, gdal, ogr
#from shapely.wkt import loads
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

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
            self.format = r'$\mathdefault{%s}$' % self.format

def map_vector_polygon_data(sFilename_in,
                            iFlag_color_in = None,
                            iFlag_colorbar_in = None,
                            sVariable_in = None,
                            sFilename_output_in=None,
                            iFlag_scientific_notation_colorbar_in=None,
                            iFont_size_in = None,
                            sColormap_in = None,
                            sTitle_in = None,
                            iDPI_in = None,
                            dMissing_value_in=None,
                            dData_max_in = None,
                            dData_min_in = None,
                            sExtend_in =None,
                            sFont_in = None,
                            sUnit_in=None,
                            aLegend_in = None,
                            aExtent_in = None,
                            pProjection_map_in=None):
    """
    plot vector data on a map
    currently only support geojson and shapefile
    by default, the program will plot all the polygons in the file
    in the furture, the program will support to plot only a subset of polygons

    Args:
        iFiletype_in (_type_): _description_
        sFilename_in (_type_): _description_
        sFilename_output_in (_type_): _description_
        iFlag_scientific_notation_colorbar_in (_type_, optional): _description_. Defaults to None.
        sColormap_in (_type_, optional): _description_. Defaults to None.
        sTitle_in (_type_, optional): _description_. Defaults to None.
        iDPI_in (_type_, optional): _description_. Defaults to None.
        dMissing_value_in (_type_, optional): _description_. Defaults to None.
        dData_max_in (_type_, optional): _description_. Defaults to None.
        dData_min_in (_type_, optional): _description_. Defaults to None.
        sExtend_in (_type_, optional): _description_. Defaults to None.
        sUnit_in (_type_, optional): _description_. Defaults to None.
        aLegend_in (_type_, optional): _description_. Defaults to None.
    """

    pDriver = ogr.GetDriverByName('GeoJSON')

    #check if the file exists
    if os.path.isfile(sFilename_in):
        pass
    else:
        print('The file does not exist: ', sFilename_in)
        return

    pDataset = pDriver.Open(sFilename_in, gdal.GA_ReadOnly)
    pLayer = pDataset.GetLayer(0)


    if iDPI_in is not None:
        iDPI = iDPI_in
    else:
        iDPI = 300

    if iFlag_color_in is not None:
        iFlag_color = iFlag_color_in
    else:
        iFlag_color = 0

    if iFlag_colorbar_in is not None:
        iFlag_colorbar = iFlag_colorbar_in
    else:
        iFlag_colorbar = 0

    if iFlag_colorbar == 1: #if colorbar is used, then the color is used
        iFlag_color = 1
    
    if iFont_size_in is not None:
        iFont_size = iFont_size_in
    else:
        iFont_size = 12

    if dMissing_value_in is not None:
        dMissing_value = dMissing_value_in
    else:
        dMissing_value = -9999

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

    if iFlag_scientific_notation_colorbar_in is not None:
        iFlag_scientific_notation_colorbar = iFlag_scientific_notation_colorbar_in
    else:
        iFlag_scientific_notation_colorbar = 0

    if sColormap_in is not None:
        sColormap = sColormap_in
    else:
        sColormap =  'rainbow'

    if sTitle_in is not None:
        sTitle = sTitle_in.capitalize()
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

    if sVariable_in is not None:
        sVariable = sVariable_in
    else:
        sVariable =  'id'

    cmap = cm.get_cmap(sColormap)

    fig = plt.figure( dpi = iDPI  )
    iSize_x= 8
    iSize_y= 8
    fig.set_figwidth( iSize_x )
    fig.set_figheight( iSize_y )

    pLayer = pDataset.GetLayer(0)
    pSrs = osr.SpatialReference()
    pSrs.ImportFromEPSG(4326)    # WGS84 lat/lon
    dLat_min = 90
    dLat_max = -90
    dLon_min = 180
    dLon_max = -180
    
    for pFeature in pLayer:
        pGeometry_in = pFeature.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()

        if sGeometry_type =='POLYGON':
            #dummy0 = loads( pGeometry_in.ExportToWkt() )            
            #aCoords_gcs = dummy0.exterior.coords
            #aCoords_gcs= np.array(aCoords_gcs)
            aCoords_gcs = get_geometry_coords(pGeometry_in)
            dLon_max = np.max( [dLon_max, np.max(aCoords_gcs[:,0])] )
            dLon_min = np.min( [dLon_min, np.min(aCoords_gcs[:,0])] )
            dLat_max = np.max( [dLat_max, np.max(aCoords_gcs[:,1])] )
            dLat_min = np.min( [dLat_min, np.min(aCoords_gcs[:,1])] )

    if iFlag_color == 1:
        aValue = list() #this one need to be fixed
        for pFeature in pLayer:    
            pGeometry_in = pFeature.GetGeometryRef()
            sGeometry_type = pGeometry_in.GetGeometryName()
            if sGeometry_type =='POLYGON':                
                dummy = float(pFeature.GetField(sVariable))
                aValue.append(dummy)

        aValue = np.array(aValue)
        if iFlag_data_min == 1  and iFlag_data_max ==1: #both are provided
            aValue = np.clip(aValue, dData_min, dData_max)
            dValue_max = dData_max # np.max(aValue)
            dValue_min = dData_min # np.min(aValue)
        else:
            aValue = aValue[aValue != dMissing_value]
            dValue_max = np.max(aValue)
            dValue_min = np.min(aValue)
            pass

        if dValue_max == dValue_min:
            iFlag_same_value = 1
            return
        else:
            iFlag_same_value = 0
            pass


    cmap = cm.get_cmap(sColormap)
    if pProjection_map_in is not None:
        pProjection_map = pProjection_map_in
    else:
        pProjection_map = ccrs.Orthographic(central_longitude =  0.50*(dLon_max+dLon_min),  central_latitude = 0.50*(dLat_max+dLat_min), globe=None)

    ax = fig.add_axes([0.08, 0.1, 0.62, 0.7], projection=pProjection_map )
    ax.set_global()
    for pFeature in pLayer:
        pGeometry_in = pFeature.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()
        #get attribute
        if iFlag_color ==1:
            dValue = float(pFeature.GetField(sVariable))
            if dValue != dMissing_value:
                if dValue > dValue_max:
                    dValue = dValue_max
                if dValue < dValue_min:
                    dValue = dValue_min

            iColor_index = int( (dValue - dValue_min) / (dValue_max - dValue_min) * 255 )
            #pick color from colormap
            sColor = cmap(iColor_index)
        else:
            sColor='black'

        if iFlag_color ==1:
            if sGeometry_type =='POLYGON':                
                #old way that has a dependency on shapely
                #dummy0 = loads( pGeometry_in.ExportToWkt() )
                #aCoords_gcs = dummy0.exterior.coords
                #aCoords_gcs= np.array(aCoords_gcs)
                #new method that use gdal directly
                aCoords_gcs = get_geometry_coords(pGeometry_in)
                polygon = mpatches.Polygon(aCoords_gcs[:,0:2], closed=True, linewidth=0.25,
                                           alpha=0.8, edgecolor = sColor,
                                           facecolor = sColor,
                                           transform=ccrs.Geodetic() )

                ax.add_patch(polygon)
            else:
                pass
        else:
            if sGeometry_type =='POLYGON':
                #dummy0 = loads( pGeometry_in.ExportToWkt() )
                #aCoords_gcs = dummy0.exterior.coords
                #aCoords_gcs= np.array(aCoords_gcs)
                aCoords_gcs = get_geometry_coords(pGeometry_in)
                polygon = mpatches.Polygon(aCoords_gcs[:,0:2],
                                           closed=True,
                                           linewidth=0.25,
                                           alpha=0.8, edgecolor = sColor,
                                           facecolor = 'none',
                                           transform= ccrs.Geodetic() )

                ax.add_patch(polygon)

    if aExtent_in is None:
        marginx  = (dLon_max - dLon_min) / 50
        marginy  = (dLat_max - dLat_min) / 50
        aExtent = [dLon_min - marginx , dLon_max + marginx , dLat_min -marginy , dLat_max + marginy]
    else:
        aExtent = aExtent_in

    ax.set_extent( aExtent )
    ax.coastlines(color='black', linewidth=1)
    ax.set_title(sTitle)
    if aLegend_in is not None:
        nlegend = len(aLegend_in)
        dLocation0 = 0.96
        for i in range(nlegend):
            sText = aLegend_in[i]
            #dLocation = 0.06 + i * 0.04
            dLocation = dLocation0 - i * 0.06
            ax.text(0.03, dLocation, sText, 
                    verticalalignment='top', horizontalalignment='left',
                    transform=ax.transAxes, 
                    color='black', fontsize=iFont_size)

            pass

    if iFlag_colorbar ==1:
        ax_cb= fig.add_axes([0.75, 0.15, 0.02, 0.6])

        if iFlag_scientific_notation_colorbar==1:
            formatter = OOMFormatter(fformat= "%1.1e")
            cb = mpl.colorbar.ColorbarBase(ax_cb, orientation='vertical',
                                           cmap=cmap,
                                           norm=mpl.colors.Normalize(dValue_min, dValue_max),  # vmax and vmin
                                           extend=sExtend, format=formatter)
        else:
            formatter = OOMFormatter(fformat= "%1.2f")
            cb = mpl.colorbar.ColorbarBase(ax_cb, orientation='vertical',
                                           cmap=cmap,
                                           norm=mpl.colors.Normalize(dValue_min, dValue_max),  # vmax and vmin
                                           extend=sExtend, format=formatter)

        cb.ax.get_yaxis().set_ticks_position('right')
        cb.ax.get_yaxis().labelpad = 5
        cb.ax.set_ylabel(sUnit, rotation=90)
        cb.ax.get_yaxis().set_label_position('left')
        cb.ax.tick_params(labelsize=6)

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=0.5, linestyle='--')
    
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    gl.xlabel_style = {'size': 10, 'color': 'k', 'rotation':0, 'ha':'right'}
    gl.ylabel_style = {'size': 10, 'color': 'k', 'rotation':90,'weight': 'normal'}
    

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
                
    #clean cache
    plt.close('all')
    plt.clf()
