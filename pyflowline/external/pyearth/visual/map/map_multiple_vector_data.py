import os
import numpy as np
from osgeo import  osr, gdal, ogr
#from shapely.wkt import loads
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.path as mpath
import matplotlib.cm as cm
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
            self.format = r'$\mathdefault{%s}$' % self.format

def map_multiple_vector_data(aFiletype_in,
                             aFilename_in,
                             iFlag_colorbar_in=None,
                             aFlag_thickness_in=None,
                             aFlag_color_in=None,
                             aFlag_fill_in = None,
                             aVariable_in = None,
                             sFilename_output_in=None,
                             iFlag_scientific_notation_colorbar_in=None,
                             iFont_size_in=None,
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

    Because of the overlay effect, it is ideal to plot them in the following order: polygon->polyling->point

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

    #check vector type first
    if aFiletype_in is None:
        print('Error: please specify the vector type')
        return
    else:
        #point: 1, polyline: 2, polygon: 3      
        arr = np.array(aFiletype_in)      
        # Check if the array is descending or has the same values
        is_descending = np.all(np.diff(arr) <= 0)        
        if is_descending == True:
            pass
        else:
            print('Error: the vector type is not correct')
            return

    pDriver = ogr.GetDriverByName('GeoJSON')


    nFile = len(aFilename_in)
    if aFlag_thickness_in is None:
        aFlag_thickness= np.zeros(nFile, dtype=np.int16)
    else:
        aFlag_thickness = aFlag_thickness_in

    if aFlag_color_in is None:
        aFlag_color= np.zeros(nFile, dtype=np.int16)
    else:
        aFlag_color = aFlag_color_in

    if aFlag_fill_in is None:
        aFlag_fill= np.zeros(nFile, dtype=np.int16)
    else:
        aFlag_fill = aFlag_fill_in

    #get the extent first
    sFilename_in = aFilename_in[0]

    pDataset = pDriver.Open(sFilename_in, gdal.GA_ReadOnly)
    pLayer = pDataset.GetLayer(0)

    if iFlag_colorbar_in is not None:
        iFlag_colorbar = iFlag_colorbar_in
    else:
        iFlag_colorbar = 0
    if iDPI_in is not None:
        iDPI = iDPI_in
    else:
        iDPI = 300

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
        sTitle = sTitle_in
        iFlag_title =1
    else:
        iFlag_title=0
        sTitle =  ''
    
    if iFont_size_in is not None:
        iFont_size = iFont_size_in
    else:
        iFont_size = 12

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
    fig = plt.figure( dpi = iDPI  )
    iSize_x= 8
    iSize_y= 8
    fig.set_figwidth( iSize_x )
    fig.set_figheight( iSize_y )


    #we require that the first polygon file defines the extent

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
        else:
            if sGeometry_type =='LINESTRING':
                aCoords_gcs = get_geometry_coords(pGeometry_in)
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

    ax = fig.add_axes([0.08, 0.1, 0.62, 0.7], projection=pProjection_map )


    #====================================
    #should we allow more than one scale for one variable?
    aValue_all = list()
    for i in range(nFile):
        aValue = list()
        sFilename = aFilename_in[i]

        iFlag_thickness = aFlag_thickness[i]
        iFlag_color = aFlag_color[i]
        if iFlag_thickness ==1 :
            sVariable = aVariable_in[i]
            pDataset = pDriver.Open(sFilename, gdal.GA_ReadOnly)
            pLayer = pDataset.GetLayer(0)

            for pFeature in pLayer:
                pGeometry_in = pFeature.GetGeometryRef()
                sGeometry_type = pGeometry_in.GetGeometryName()
                dValue = float(pFeature.GetField(sVariable))
                aValue.append(dValue)
        else:
            if iFlag_color  == 1:
                sVariable = aVariable_in[i]
                pDataset = pDriver.Open(sFilename, gdal.GA_ReadOnly)
                pLayer = pDataset.GetLayer(0)

                for pFeature in pLayer:
                    pGeometry_in = pFeature.GetGeometryRef()
                    sGeometry_type = pGeometry_in.GetGeometryName()
                    dValue = float(pFeature.GetField(sVariable))
                    aValue.append(dValue)


        aValue_all.append(aValue)

    iThickness_max = 2.5
    iThickness_min = 0.3

    for i in range(nFile):
        sFilename = aFilename_in[i]
        iFlag_thickness = aFlag_thickness[i]
        iFlag_color = aFlag_color[i]
        iFlag_fill = aFlag_fill[i]
        pDataset = pDriver.Open(sFilename, gdal.GA_ReadOnly)
        pLayer = pDataset.GetLayer(0)

        nColor = pLayer.GetFeatureCount()
        aColor = cm.rainbow(np.linspace(0, 1, nColor))
        lID = 0
        for pFeature in pLayer:
            pGeometry_in = pFeature.GetGeometryRef()
            sGeometry_type = pGeometry_in.GetGeometryName()

            if iFlag_thickness ==1 :
                sVariable = aVariable_in[i]
                aValue = np.array(aValue_all[i])
                if iFlag_data_min == 1  and iFlag_data_max ==1: #both are provided
                    aValue = np.clip(aValue, dData_min, dData_max)
                    dValue_max = dData_max
                    dValue_min = dData_min
                else:
                    aValue = aValue[aValue != dMissing_value]
                    dValue_max = np.max(aValue)
                    dValue_min = np.min(aValue)
                    pass
                dValue = float(pFeature.GetField(sVariable))

                if dValue != dMissing_value:
                    if dValue > dValue_max:
                        dValue = dValue_max

                    if dValue < dValue_min:
                        dValue = dValue_min

                iThickness = remap( dValue, dValue_min, dValue_max, iThickness_min, iThickness_max )
            else:
                iThickness = 0.25

            if iFlag_color ==1:
                if nColor < 10:
                    sColor = aColor[lID]
                else:
                    sVariable = aVariable_in[i]
                    aValue = np.array(aValue_all[i])
                    if iFlag_data_min == 1  and iFlag_data_max ==1: #both are provided
                        aValue = np.clip(aValue, dData_min, dData_max)
                        dValue_max = dData_max
                        dValue_min = dData_min
                    else:
                        aValue = aValue[aValue != dMissing_value]
                        dValue_max = np.max(aValue)
                        dValue_min = np.min(aValue)
                        pass
                    dValue = float(pFeature.GetField(sVariable))

                    if dValue != dMissing_value:
                        if dValue > dValue_max:
                            dValue = dValue_max

                        if dValue < dValue_min:
                            dValue = dValue_min

                    iColor_index = int( (dValue - dValue_min) / (dValue_max - dValue_min) * 255 )
                    sColor = cmap(iColor_index)
            else:
                sColor = 'black'

                #pick color from colormap

            if sGeometry_type =='POINT':
                #dummy0 = loads( pGeometry_in.ExportToWkt() )
                #aCoords_gcs = dummy0.coords
                #aCoords_gcs= np.array(aCoords_gcs)
                aCoords_gcs = get_geometry_coords(pGeometry_in)
                aCoords_gcs = aCoords_gcs[:,0:2]
                ax.plot(aCoords_gcs[0], aCoords_gcs[1], 'o', color= sColor, markersize=2, transform=ccrs.Geodetic())                
            else:
                if sGeometry_type =='LINESTRING':
                    #dummy0 = loads( pGeometry_in.ExportToWkt() )
                    #aCoords_gcs = dummy0.coords
                    #aCoords_gcs= np.array(aCoords_gcs)
                    aCoords_gcs = get_geometry_coords(pGeometry_in)
                    aCoords_gcs = aCoords_gcs[:,0:2]
                    nvertex = len(aCoords_gcs)
                    codes = np.full(nvertex, mpath.Path.LINETO, dtype=int )
                    codes[0] = mpath.Path.MOVETO
                    path = mpath.Path(aCoords_gcs, codes)
                    x, y = zip(*path.vertices)
                    line, = ax.plot(x, y, color= sColor, linewidth=iThickness, transform=ccrs.Geodetic())
                else:
                    if sGeometry_type == 'POLYGON': 
                        #dummy0 = loads( pGeometry_in.ExportToWkt() )
                        #aCoords_gcs = dummy0.exterior.coords
                        #aCoords_gcs = np.array(aCoords_gcs)
                        aCoords_gcs = get_geometry_coords(pGeometry_in)
                        if iFlag_fill == 1:
                            polygon = mpatches.Polygon(aCoords_gcs[:,0:2], closed=True, linewidth=0.25, 
                                           alpha=0.8, edgecolor = sColor,facecolor= sColor, 
                                           transform=ccrs.Geodetic() )
                        else:
                            polygon = mpatches.Polygon(aCoords_gcs[:,0:2], closed=True, linewidth=0.25, 
                                           alpha=0.8, edgecolor = sColor,facecolor='none', 
                                           transform=ccrs.Geodetic() )
                        ax.add_patch(polygon)
                    else:
                        pass

            lID = lID + 1

    if aExtent_in is None:
        marginx  = (dLon_max - dLon_min) / 50
        marginy  = (dLat_max - dLat_min) / 50
        aExtent = [dLon_min - marginx , dLon_max + marginx , dLat_min -marginy , dLat_max + marginy]
    else:
        aExtent = aExtent_in

    ax.set_global()
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
            formatter = OOMFormatter(fformat= "%1.1f")
            cb = mpl.colorbar.ColorbarBase(ax_cb, orientation='vertical',
                                           cmap=cmap,
                                           norm=mpl.colors.Normalize(dValue_min, dValue_max),  # vmax and vmin
                                           extend=sExtend, format=formatter)

        cb.ax.get_yaxis().set_ticks_position('right')
        cb.ax.get_yaxis().labelpad = 10
        cb.ax.set_ylabel(sUnit, rotation=270)
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
