import os
import json
#dependency packages
import numpy as np
from osgeo import  osr, gdal, ogr
from shapely.wkt import loads
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches
import matplotlib.cm as cm

from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.io.shapereader as shpreader
from cartopy.feature import ShapelyFeature
import cartopy.crs as ccrs

import requests
#import utm

from pyflowline.algorithms.auxiliary.text_reader_string import text_reader_string
from pyflowline.algorithms.auxiliary.gdal_functions import Google_MetersPerPixel
from pyflowline.algorithms.auxiliary.gdal_functions import reproject_coordinates, reproject_coordinates_batch
from pyflowline.formats.convert_flowline_to_geojson import convert_shapefile_to_geojson_swat
from pyflowline.algorithms.auxiliary.gdal_functions import gdal_read_geotiff_file, retrieve_geotiff_metadata

#the orthographic projection should be used on a sphere, this might be provided as an input
pProjection_map = ccrs.Orthographic(central_longitude =  0.50*(-149.5+(-146.5)), \
        central_latitude = 0.50*(68.1+70.35), globe=None)
        

def _plot(self, sFilename_in,iFlag_title=None, sVariable_in=None, aExtent_in = None, pProjection_map_in = None):
    if sVariable_in == 'mesh':
        self._plot_mesh(sFilename_in, aExtent_in= aExtent_in, pProjection_map_in= pProjection_map_in)
    else:            
        if sVariable_in == 'overlap':
            self._plot_mesh_with_flowline( sFilename_in, iFlag_title= iFlag_title, aExtent_in= aExtent_in, pProjection_map_in= pProjection_map_in)
        else:            
            
            for pBasin in self.aBasin:            
                pBasin._basin_plot(self.iCase_index, self.iMesh_type, self.sMesh_type,sFilename_in,\
                    iFlag_title= iFlag_title, \
                    sVariable_in= sVariable_in, \
                        pProjection_map_in= pProjection_map_in)
            pass
    
    return
    
def _plot_study_area(self, sFilename_boundary_in = None, sFilename_slope_in = None, sFilename_nhd_in = None):
    sWorkspace_output_case = self.sWorkspace_output
    fig = plt.figure(dpi=300)
    fig.set_figwidth( 4 )
    fig.set_figheight( 4 )
    #plot dem first        
    sFilename_dem= self.sFilename_dem
     
    dummy = gdal_read_geotiff_file(sFilename_dem)
    dResolution_x = dummy[1]
    pProjection = dummy[8]
    pSpatial_reference_source = dummy[9]
    dem_proj =ccrs.AlbersEqualAre(central_longitude=pSpatial_reference_source.GetProjPar('longitude_of_center'), \
        central_latitude=pSpatial_reference_source.GetProjPar('latitude_of_center'),\
            standard_parallels=(pSpatial_reference_source.GetProjPar('standard_parallel_1'),\
                pSpatial_reference_source.GetProjPar('standard_parallel_2')))
    ax_dem = fig.add_axes([0.3, 0.15, 0.6, 0.55] , projection=dem_proj )
    ax_dem.set_xmargin(0.05)
    ax_dem.set_ymargin(0.10)          
    
    missing_value = dummy[6]
    
    dOriginX=dummy[2]
    dOriginY=dummy[3]
    nrow=dummy[4]
    ncolumn=dummy[5]
    dLon_min = dOriginX
    dLon_max = dOriginX + ncolumn * dResolution_x
    dLat_max = dOriginY
    dLat_min = dOriginY - nrow * dResolution_x
   
    pSpatial_reference_target = osr.SpatialReference()  
    pSpatial_reference_target.ImportFromEPSG(4326)
    aLon= list()
    aLat=list()
    aLon.append(dLon_min)
    aLon.append(dLon_max)
    aLat.append(dLat_min)
    aLat.append(dLat_max)
    aLon, aLat = reproject_coordinates_batch(aLon, aLat,pSpatial_reference_source, pSpatial_reference_target)
    dLongitude_center = np.mean(aLon)
    dLatitude_center = np.mean(aLat)
    aImage_extent =  [dLon_min- dResolution_x ,dLon_max + dResolution_x,dLat_min -dResolution_x,  dLat_max+dResolution_x]
    aImage_in = dummy[0]
    aImage_in[np.where(aImage_in == missing_value)] = np.nan
    # 
    demplot = ax_dem.imshow(aImage_in, origin='upper', extent=aImage_extent,cmap=cm.terrain,transform=dem_proj) 
    # change all spines
    for axis in ['top','bottom','left','right']:
        ax_dem.spines[axis].set_linewidth(0.5)
    # increase tick width
    ax_dem.tick_params(width=0.5)
    # unused code
    # u = utm.from_latlon(dLatitude_center, dLongitude_center)
    gl = ax_dem.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.3, linestyle='--')

    gl.left_labels=False
    gl.top_labels=False        
   
     # declare text size
    XTEXT_SIZE = 8
    YTEXT_SIZE = 8
    # to facilitate text rotation at bottom edge, ...
    # text justification: 'ha':'right' is used to avoid clashing with map'sboundary
    # default of 'ha' is center, often causes trouble when text rotation isnot zero
    gl.xlabel_style = {'size': XTEXT_SIZE, 'color': 'k', 'rotation':0,'ha':'right'}
    gl.ylabel_style = {'size': YTEXT_SIZE, 'color': 'k', 'rotation':90, 'weight': 'normal'}

    
    if sFilename_nhd_in is None:
        sFilename_nhd = '/qfs/people/liao313/data/hexwatershed/susquehannavector/hydrology/nhd_proj.shp'
    else:
        sFilename_nhd = sFilename_nhd_in
    
    reader = shpreader.Reader(sFilename_nhd)
    shape_feature = ShapelyFeature(reader.geometries(),
                            dem_proj, facecolor='none')
    ax_dem.add_feature(shape_feature)
    pSpatial_reference_source = osr.SpatialReference()  
    pSpatial_reference_source.ImportFromEPSG(4326)
    pSpatial_reference_target = osr.SpatialReference()  
    pSpatial_reference_target.ImportFromWkt(dem_proj.to_wkt())
    for pBasin in self.aBasin:
        sFilename_dam = pBasin.sFilename_dam
        aData_dam = text_reader_string(sFilename_dam, iSkipline_in =1, cDelimiter_in=',' )
        ndam = len(aData_dam)
        for idam in range(ndam):
            dLon = float(aData_dam[idam, 1])
            dLat = float(aData_dam[idam, 0])
            x,y = reproject_coordinates(dLon, dLat, pSpatial_reference_source, pSpatial_reference_target )
            ax_dem.plot(x,y,marker = 'x', color='red', markersize=4)
            pass

    ax_cb= fig.add_axes([0.2, 0.2, 0.02, 0.5])    
    cb = plt.colorbar(demplot, cax = ax_cb, extend = 'both')
    cb.ax.get_yaxis().set_ticks_position('left')
    cb.ax.get_yaxis().labelpad = 10
    cb.ax.set_ylabel('Unit: meter', rotation=270)
    cb.ax.tick_params(labelsize=6) 
    #google earth   
    ge_proj = ccrs.Mercator()#central_longitude=dLongitude_center)
    
    pSpatial_reference_source = osr.SpatialReference()  
    pSpatial_reference_source.ImportFromEPSG(4326)
    pSpatial_reference_target = osr.SpatialReference()  
    
    iZoom = 7       
    pSpatial_reference_target.ImportFromWkt(ge_proj.to_wkt())
    uv_xcenter, uv_ycenter = reproject_coordinates(dLongitude_center,dLatitude_center, pSpatial_reference_source, pSpatial_reference_target)
    xsize_ge = 600
    ysize_ge = 600
    scale = Google_MetersPerPixel(iZoom)
    xrange = (uv_xcenter - (xsize_ge/2.0*scale), uv_xcenter + (xsize_ge/20*scale))
    yrange = (uv_ycenter - (ysize_ge/2.0*scale), uv_ycenter + (ysize_ge/20*scale))
    aImage_extent = [xrange[0], xrange[1],yrange[0], yrange[1] ]
    ax_ge = fig.add_axes([0.1, 0.75, 0.2, 0.2] , projection=ge_proj )
     
    sLongitude_center ="{:0f}".format(dLongitude_center)
    sLatitude_center= "{:0f}".format(dLatitude_center)
    sZoom="{:0d}".format(iZoom)
    
    sResolution = "{:0d}".format(xsize_ge) + 'x' + "{:0d}".format(ysize_ge)
    sMap_type="hybrid"
    sGoogleMap = "http://maps.googleapis.com/maps/api/staticmap?" + \
        "center=" + sLatitude_center + ',' + sLongitude_center + \
        "&zoom=" + sZoom + "&size=" + sResolution + \
        "&maptype="+sMap_type+"&sensor=false&format=png32" + "key=AIzaSyCT0NPGHJqRFysOYELOWBdqov6AbphgaFY"
    r = requests.get(sGoogleMap)  
    # wb mode is stand for write binary mode
    f = open('google_map.png', 'wb')
    # r.content gives content,
    # in this case gives image
    f.write(r.content)
    # close method of file object
    # save and close the file
    f.close()
    
    img = mpimg.imread('google_map.png')
    ax_ge.imshow(img,extent=aImage_extent,       transform=ge_proj)
    gl_ge = ax_ge.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.3, linestyle='--')
    #gl_ge.xlocator = mticker.FixedLocator(aLon_range)
    #gl_ge.ylocator = mticker.FixedLocator(aLat_range)
    gl_ge.xformatter = LongitudeFormatter()
    gl_ge.yformatter = LatitudeFormatter()
    gl_ge.right_labels=False
    gl_ge.bottom_labels=False        
   
    # declare text size
    XTEXT_SIZE = 2.5
    YTEXT_SIZE = 2.5
    # to facilitate text rotation at bottom edge, ...
    # text justification: 'ha':'right' is used to avoid clashing with map'sboundary
    # default of 'ha' is center, often causes trouble when text rotation isnot zero
    gl_ge.xlabel_style = {'size': XTEXT_SIZE, 'color': 'k', 'rotation':0,'ha':'right'}
    gl_ge.ylabel_style = {'size': YTEXT_SIZE, 'color': 'k', 'rotation':90,'weight': 'normal'}
    #add boundary
    if sFilename_boundary_in is None:
        sFilename_boundary = '/qfs/people/liao313/data/hexwatershed/susquehannavector/hydrology/boundary_ge.shp'
    else:
        sFilename_boundary = sFilename_boundary_in

    reader = shpreader.Reader(sFilename_boundary)
    shape_feature = ShapelyFeature(reader.geometries(),
                            ge_proj, facecolor='none')
    ax_ge.add_feature(shape_feature)
    #histogram
    ax_histo = fig.add_axes([0.4, 0.8, 0.4, 0.15]  )
    if sFilename_slope_in is None:
        sFilename_slope = '/qfs/people/liao313/data/hexwatershed/susquehannaraster/dem/slope_ext.tif'
    else:
        sFilename_slope= sFilename_slope_in

    dummy = gdal_read_geotiff_file(sFilename_slope)
    aSlope= dummy[0]
    missing_value = dummy[6]
    aSlope = aSlope[np.where(aSlope != missing_value)]
    dMax_x=60
    dMin_x=0
    dSpace_x=4
    ax_histo.hist(aSlope,  int(  (dMax_x-dMin_x)/ dSpace_x), color ="skyblue", ec="skyblue")  
    sLabel_x= 'Slope (degree)'
    sLabel_y ='Frequency'
    ax_histo.set_xlabel(sLabel_x,fontsize=4 )
    ax_histo.set_ylabel(sLabel_y,fontsize=4 )       
    ax_histo.set_xlim( dMin_x, dMax_x )
    formatter = mticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)        
    ax_histo.yaxis.set_major_formatter(formatter)
    ax_histo.tick_params(axis='x', labelsize=5 )
    ax_histo.tick_params(axis='y', labelsize=5 )
    ax_histo.yaxis.offsetText.set_fontsize(5)
    sText =  'Outlet'
    ax_dem.text(0.85, 0.05, sText, \
            verticalalignment='top', horizontalalignment='left',\
            transform=ax_dem.transAxes, \
            color='black', fontsize=7)
    sFilename  = 'study_area.png'
    sFilename_out = os.path.join(sWorkspace_output_case, sFilename)
    plt.savefig(sFilename_out, bbox_inches='tight')
    
    return

def _plot_mesh(self, sFilename_in, aExtent_in=None, pProjection_map_in = None):
    sWorkspace_output_case = self.sWorkspace_output
    sFilename_json  =  self.sFilename_mesh
    sMesh_type = self.sMesh_type
    fig = plt.figure(dpi=300)
    fig.set_figwidth( 4 )
    fig.set_figheight( 4 )
    ax = fig.add_axes([0.1, 0.15, 0.75, 0.7] , projection=pProjection_map )
    ax.set_global()
    pDriver = ogr.GetDriverByName('GeoJSON')
    pDataset = pDriver.Open(sFilename_json, gdal.GA_ReadOnly)
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
        lID =0 
        if sGeometry_type =='POLYGON':
            dummy0 = loads( pGeometry_in.ExportToWkt() )
            aCoords_gcs = dummy0.exterior.coords
            aCoords_gcs= np.array(aCoords_gcs)
            nvertex = len(aCoords_gcs)
            for i in range(nvertex):
                dLon = aCoords_gcs[i][0]
                dLat = aCoords_gcs[i][1]
                if dLon > dLon_max:
                    dLon_max = dLon
                if dLon < dLon_min:
                    dLon_min = dLon
                if dLat > dLat_max:
                    dLat_max = dLat
                if dLat < dLat_min:
                    dLat_min = dLat
            polygon = mpatches.Polygon(aCoords_gcs[:,0:2], closed=True, linewidth=0.25, \
                alpha=0.8, edgecolor = 'black',facecolor='none', \
                    transform=ccrs.PlateCarree() )
            ax.add_patch(polygon)                   
    
    if aExtent_in is None:
        marginx  = (dLon_max - dLon_min) / 20
        marginy  = (dLat_max - dLat_min) / 20
        aExtent = [dLon_min - marginx , dLon_max + marginx , dLat_min -marginy , dLat_max + marginy]
    else:
        aExtent = aExtent_in
    
    ax.set_extent( aExtent )  

    #add mesh info
    ax.set_title( sMesh_type.title() +  ' mesh')
    
    ax.coastlines()#resolution='110m')
    ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.3, linestyle='--')
    sDirname = os.path.dirname(sFilename_json)
    #sFilename  = Path(sFilename_json).stem + '.png'
    sFilename_out = os.path.join(sDirname, sFilename_in)
    plt.savefig(sFilename_out, bbox_inches='tight')
    pDataset = pLayer = pFeature  = None   
      
    return

def _plot_mesh_with_flowline(self, sFilename_in, iFlag_title=None, aExtent_in=None, pProjection_map_in = None):
    sWorkspace_output_case = self.sWorkspace_output
    sFilename_mesh  =  self.sFilename_mesh
    sMesh_type = self.sMesh_type
    fig = plt.figure( dpi=300 )
    fig.set_figwidth( 4 )
    fig.set_figheight( 4 )
    ax = fig.add_axes([0.1, 0.15, 0.75, 0.7] , projection=pProjection_map )
    ax.set_global()
    pDriver = ogr.GetDriverByName('GeoJSON')
    pDataset = pDriver.Open(sFilename_mesh, gdal.GA_ReadOnly)
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
            dummy0 = loads( pGeometry_in.ExportToWkt() )
            aCoords_gcs = dummy0.exterior.coords
            aCoords_gcs= np.array(aCoords_gcs)
            nvertex = len(aCoords_gcs)
            for i in range(nvertex):
                dLon = aCoords_gcs[i][0]
                dLat = aCoords_gcs[i][1]
                if dLon > dLon_max:
                    dLon_max = dLon
                if dLon < dLon_min:
                    dLon_min = dLon
                if dLat > dLat_max:
                    dLat_max = dLat
                if dLat < dLat_min:
                    dLat_min = dLat
            polygon = mpatches.Polygon(aCoords_gcs[:,0:2], closed=True,  linewidth=0.25, \
                alpha=0.8, edgecolor = 'black',facecolor='none', \
                    transform=ccrs.PlateCarree() )
            ax.add_patch(polygon)    
                           
    #draw base flowline first with black color
    lID = 0 
    for pBasin in self.aBasin:
        sWorkspace_output_basin=  pBasin.sWorkspace_output_basin                        
        sFilename = pBasin.sFilename_flowline_simplified
        sFilename_out = os.path.join(sWorkspace_output_basin, sFilename)
        sFilename_json = os.path.join(sWorkspace_output_basin,sFilename_out)
        pDriver = ogr.GetDriverByName('GeoJSON')
        pDataset = pDriver.Open(sFilename_json, gdal.GA_ReadOnly)
        pLayer = pDataset.GetLayer(0)
        n_colors = pLayer.GetFeatureCount()     
        for pFeature in pLayer:
            pGeometry_in = pFeature.GetGeometryRef()
            sGeometry_type = pGeometry_in.GetGeometryName()
            if sGeometry_type =='LINESTRING':
                dummy0 = loads( pGeometry_in.ExportToWkt() )
                aCoords_gcs = dummy0.coords
                aCoords_gcs= np.array(aCoords_gcs)
                nvertex = len(aCoords_gcs)   
                codes = np.full(nvertex, mpath.Path.LINETO, dtype=int )
                codes[0] = mpath.Path.MOVETO
                path = mpath.Path(aCoords_gcs, codes)            
                x, y = zip(*path.vertices)
                line, = ax.plot(x, y, color= 'black', linewidth=0.5)
                lID = lID + 1
            pass
        pass
    #plot new flowline now
    lID = 0 
    for pBasin in self.aBasin:
        sWorkspace_output_basin=  pBasin.sWorkspace_output_basin            
        sFilename_out = pBasin.sFilename_flowline_conceptual           
        sFilename_json = os.path.join(sWorkspace_output_basin,sFilename_out)
        pDriver = ogr.GetDriverByName('GeoJSON')
        pDataset = pDriver.Open(sFilename_json, gdal.GA_ReadOnly)
        pLayer = pDataset.GetLayer(0)
        n_colors = pLayer.GetFeatureCount()
    
        colours = cm.rainbow(np.linspace(0, 1, n_colors))
        for pFeature in pLayer:
            pGeometry_in = pFeature.GetGeometryRef()
            sGeometry_type = pGeometry_in.GetGeometryName()
            if sGeometry_type =='LINESTRING':
                dummy0 = loads( pGeometry_in.ExportToWkt() )
                aCoords_gcs = dummy0.coords
                aCoords_gcs= np.array(aCoords_gcs)
                nvertex = len(aCoords_gcs)   
                codes = np.full(nvertex, mpath.Path.LINETO, dtype=int )
                codes[0] = mpath.Path.MOVETO
                path = mpath.Path(aCoords_gcs, codes)            
                x, y = zip(*path.vertices)
                line, = ax.plot(x, y, color= colours[lID], linewidth=1)
                lID = lID + 1
            pass
        pass
    
    #dam
    pSpatial_reference_source = osr.SpatialReference()  
    pSpatial_reference_source.ImportFromEPSG(4326)
    pSpatial_reference_target = osr.SpatialReference() 
    for pBasin in self.aBasin:
        sFilename_dam = pBasin.sFilename_dam
        iFlag_dam= pBasin.iFlag_dam
        if iFlag_dam==1:
            aData_dam = text_reader_string(sFilename_dam, iSkipline_in =1, cDelimiter_in=',' )
            ndam = len(aData_dam)
            for idam in range(ndam):
                dLon = float(aData_dam[idam, 1])
                dLat = float(aData_dam[idam, 0])                
                ax.plot(dLon,dLat,marker = 'x', color='red', markersize=4)
                pass    
    sDirname = os.path.dirname(sFilename_mesh)
    if aExtent_in is None:
        marginx  = (dLon_max - dLon_min) / 20
        marginy  = (dLat_max - dLat_min) / 20
        aExtent = [dLon_min - marginx , dLon_max + marginx , dLat_min -marginy , dLat_max + marginy]           
    else:
        aExtent = aExtent_in       
    
    sTitle = 'Conceptual flowline'                     
    
    ax.set_extent( aExtent )        
    ax.coastlines()
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=0.2, color='gray', alpha=0.3, linestyle='--')
    gl.xlabel_style = {'size': 8, 'color': 'k', 'rotation':0, 'ha':'right'}
    gl.ylabel_style = {'size': 8, 'color': 'k', 'rotation':90,'weight':'normal'}
    gl.xlocator = mticker.MaxNLocator(5)
    gl.ylocator = mticker.MaxNLocator(5)

    if iFlag_title is None:
        ax.set_title( sTitle )
    else:
        if iFlag_title==1:
            ax.set_title( sTitle )
        else:
            pass
    iCase = self.iCase_index
    if self.sMesh_type == 'hexagon':
        iCase = 10 - iCase
    if self.sMesh_type == 'square':
        iCase = 10 - iCase
    if self.sMesh_type == 'latlon':
        iCase = 10 - iCase 
    sText = 'Case: ' + "{:0d}".format( int(iCase) ) 
    ax.text(0.05, 0.95, sText, \
    verticalalignment='top', horizontalalignment='left',\
            transform=ax.transAxes, \
            color='black', fontsize=8)

    sText = 'Mesh type: ' + sMesh_type.title()
    ax.text(0.05, 0.9, sText, \
    verticalalignment='top', horizontalalignment='left',\
            transform=ax.transAxes, \
            color='black', fontsize=8)
    
    sResolution =  'Resolution: ' + "{:0d}".format( int(self.dResolution_meter) ) + 'm'
    if self.sMesh_type != 'mpas':
        ax.text(0.05, 0.85, sResolution, \
            verticalalignment='top', horizontalalignment='left',\
            transform=ax.transAxes, \
            color='black', fontsize=8)
    
    
    sFilename_out = os.path.join(sDirname, sFilename_in)
    plt.savefig(sFilename_out, bbox_inches='tight')
    pDataset = pLayer = pFeature  = None   
    
    plt.close(fig)
    return

def _compare_with_raster_dem_method(self, sFilename_dem_flowline, sFilename_in, aExtent_in=None, pProjection_map_in = None):
    sWorkspace_output_case = self.sWorkspace_output
    sFilename_mesh  =  self.sFilename_mesh
    sMesh_type = self.sMesh_type
    fig = plt.figure( dpi=300 )
    fig.set_figwidth( 4 )
    fig.set_figheight( 4 )
    ax = fig.add_axes([0.1, 0.15, 0.75, 0.7] , projection=pProjection_map )
    ax.set_global()
    pDriver = ogr.GetDriverByName('GeoJSON')
    pDataset = pDriver.Open(sFilename_mesh, gdal.GA_ReadOnly)
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
            dummy0 = loads( pGeometry_in.ExportToWkt() )
            aCoords_gcs = dummy0.exterior.coords
            aCoords_gcs= np.array(aCoords_gcs)
            nvertex = len(aCoords_gcs)
            for i in range(nvertex):
                dLon = aCoords_gcs[i][0]
                dLat = aCoords_gcs[i][1]
                if dLon > dLon_max:
                    dLon_max = dLon
                if dLon < dLon_min:
                    dLon_min = dLon
                if dLat > dLat_max:
                    dLat_max = dLat
                if dLat < dLat_min:
                    dLat_min = dLat
            polygon = mpatches.Polygon(aCoords_gcs[:,0:2], closed=True,  linewidth=0.25, \
                alpha=0.8, edgecolor = 'black',facecolor='none', \
                    transform=ccrs.PlateCarree() )
            ax.add_patch(polygon)                   
    #draw base flowline first with black color
    lID = 0 
    for pBasin in self.aBasin:
        sWorkspace_output_basin=  pBasin.sWorkspace_output_basin                        
        sFilename = pBasin.sFilename_flowline_simplified
        sFilename_out = os.path.join(sWorkspace_output_basin, sFilename)
        sFilename_json = os.path.join(sWorkspace_output_basin,sFilename_out)
        pDriver = ogr.GetDriverByName('GeoJSON')
        pDataset = pDriver.Open(sFilename_json, gdal.GA_ReadOnly)
        pLayer = pDataset.GetLayer(0)
        n_colors = pLayer.GetFeatureCount()     
        for pFeature in pLayer:
            pGeometry_in = pFeature.GetGeometryRef()
            sGeometry_type = pGeometry_in.GetGeometryName()
            if sGeometry_type =='LINESTRING':
                dummy0 = loads( pGeometry_in.ExportToWkt() )
                aCoords_gcs = dummy0.coords
                aCoords_gcs= np.array(aCoords_gcs)
                nvertex = len(aCoords_gcs)   
                codes = np.full(nvertex, mpath.Path.LINETO, dtype=int )
                codes[0] = mpath.Path.MOVETO
                path = mpath.Path(aCoords_gcs, codes)            
                x, y = zip(*path.vertices)
                line, = ax.plot(x, y, color= 'black', linewidth=0.5)
                lID = lID + 1
            pass
        pass
    #plot new flowline now
    lID = 0 
    for pBasin in self.aBasin:
        sWorkspace_output_basin=  pBasin.sWorkspace_output_basin            
        sFilename_out = pBasin.sFilename_flowline_conceptual           
        sFilename_json = os.path.join(sWorkspace_output_basin,sFilename_out)
        pDriver = ogr.GetDriverByName('GeoJSON')
        pDataset = pDriver.Open(sFilename_json, gdal.GA_ReadOnly)
        pLayer = pDataset.GetLayer(0)
        
        for pFeature in pLayer:
            pGeometry_in = pFeature.GetGeometryRef()
            sGeometry_type = pGeometry_in.GetGeometryName()
            if sGeometry_type =='LINESTRING':
                dummy0 = loads( pGeometry_in.ExportToWkt() )
                aCoords_gcs = dummy0.coords
                aCoords_gcs= np.array(aCoords_gcs)
                nvertex = len(aCoords_gcs)   
                codes = np.full(nvertex, mpath.Path.LINETO, dtype=int )
                codes[0] = mpath.Path.MOVETO
                path = mpath.Path(aCoords_gcs, codes)            
                x, y = zip(*path.vertices)
                line, = ax.plot(x, y, color= 'blue', linewidth=1)
                lID = lID + 1
            pass
        pass
    #now swat result
    #conver file first
    for pBasin in self.aBasin:
        sWorkspace_output_basin=  pBasin.sWorkspace_output_basin 
        sFilename_swat = os.path.join(sWorkspace_output_basin, 'swat.json')
        convert_shapefile_to_geojson_swat(1, sFilename_dem_flowline,sFilename_swat)
        pDriver = ogr.GetDriverByName('GeoJSON')
        pDataset = pDriver.Open(sFilename_swat, gdal.GA_ReadOnly)
        pLayer = pDataset.GetLayer(0)
       
        for pFeature in pLayer:
            pGeometry_in = pFeature.GetGeometryRef()
            sGeometry_type = pGeometry_in.GetGeometryName()
            if sGeometry_type =='LINESTRING':
                dummy0 = loads( pGeometry_in.ExportToWkt() )
                aCoords_gcs = dummy0.coords
                aCoords_gcs= np.array(aCoords_gcs)
                nvertex = len(aCoords_gcs)   
                codes = np.full(nvertex, mpath.Path.LINETO, dtype=int )
                codes[0] = mpath.Path.MOVETO
                path = mpath.Path(aCoords_gcs, codes)            
                x, y = zip(*path.vertices)
                line, = ax.plot(x, y, color= 'red', linewidth=1)
                lID = lID + 1
            pass
        pass
        #plot new
    
    sDirname = os.path.dirname(sFilename_mesh)
    if aExtent_in is None:
        marginx  = (dLon_max - dLon_min) / 20
        marginy  = (dLat_max - dLat_min) / 20
        aExtent = [dLon_min - marginx , dLon_max + marginx , dLat_min -marginy , dLat_max + marginy]     
    else:
        aExtent = aExtent_in       
   
    sTitle = 'Conceptual flowline'                 
    
    ax.set_extent( aExtent )        
    ax.coastlines()#resolution='110m')
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=0.2, color='gray', alpha=0.3, linestyle='--')
    gl.xlabel_style = {'size': 8, 'color': 'k', 'rotation':0, 'ha':'right'}
    gl.ylabel_style = {'size': 8, 'color': 'k', 'rotation':90,'weight':'normal'}
    gl.xlocator = mticker.MaxNLocator(5)
    gl.ylocator = mticker.MaxNLocator(5)
    ax.set_title( sTitle )
    sText = 'Mesh type: ' + sMesh_type.title()
    ax.text(0.05, 0.95, sText, \
    verticalalignment='top', horizontalalignment='left',\
            transform=ax.transAxes, \
            color='black', fontsize=8)
    
    sResolution =  'Resolution: ' + "{:0d}".format( int(self.dResolution_meter/1000) ) + ' km'
    if self.sMesh_type != 'mpas':
        ax.text(0.05, 0.90, sResolution, \
            verticalalignment='top', horizontalalignment='left',\
            transform=ax.transAxes, \
            color='black', fontsize=8)
    
    sText = 'Topologic relationship-based'
    ax.text(0.05, 0.85, sText, \
    verticalalignment='top', horizontalalignment='left',\
            transform=ax.transAxes, \
            color='blue', fontsize=8)
    
    sText = 'Raster DEM-based'
    ax.text(0.05, 0.80, sText, \
    verticalalignment='top', horizontalalignment='left',\
            transform=ax.transAxes, \
            color='red', fontsize=8)
    
    
    sFilename_out = os.path.join(self.sWorkspace_output, sFilename_in)
    plt.savefig(sFilename_out, bbox_inches='tight')
    pDataset = pLayer = pFeature  = None   
    plt.close(fig)
    return
  
def _basinplot(self, iCase_index, iMesh_type, sMesh_type, sFilename_in, iFlag_title=None, sVariable_in=None, aExtent_in = None, pProjection_map_in = None):
    iFlag_label = 0
    sWorkspace_output_basin = self.sWorkspace_output_basin
    if sVariable_in is not None:
        if sVariable_in == 'flowline_raw':
            sFilename_json = self.sFilename_flowline_raw
            sTitle = 'Original flowline'
        else:
            if sVariable_in == 'flowline_filter':
                sFilename_json = self.sFilename_flowline_filter
                sTitle = 'Filtered flowline'
            else:
                if sVariable_in == 'flowline_simplified':
                    sFilename_out = self.sFilename_flowline_simplified
                    sFilename_json = os.path.join(sWorkspace_output_basin, sFilename_out)
                    sTitle = 'Simplified flowline'
                    iFlag_label = 1
                else:
                    if sVariable_in == 'flowline_conceptual':
                        sFilename_out = self.sFilename_flowline_conceptual
                        sFilename_json = os.path.join(sWorkspace_output_basin, sFilename_out)
                        sTitle = 'Conceptual flowline'
                        iFlag_label = 1
                    else:
                        if sVariable_in == 'aof':
                            sFilename_out = 'area_of_difference.geojson'
                            sFilename_json = os.path.join(sWorkspace_output_basin, sFilename_out)
                            sTitle = 'Conceptual flowline'
                            iFlag_label = 1
                            self.plot_area_of_difference(iCase_index, iMesh_type, sMesh_type, sFilename_in, aExtent_in = aExtent_in)
                            return
                        else:
                            pass
                pass
    else:
        #default 
        sFilename_json = self.sFilename_flowline_conceptual
    
    #request = cimgt.OSM()
    fig = plt.figure( dpi=300)
    fig.set_figwidth( 4 )
    fig.set_figheight( 4 )
    ax = fig.add_axes([0.1, 0.15, 0.75, 0.8] , projection=pProjection_map ) #request.crs
    
    pDriver = ogr.GetDriverByName('GeoJSON')
    pDataset = pDriver.Open(sFilename_json, gdal.GA_ReadOnly)
    pLayer = pDataset.GetLayer(0)

    pSrs = osr.SpatialReference()  
    pSrs.ImportFromEPSG(4326)    # WGS84 lat/lon

    lID = 0
    dLat_min = 90
    dLat_max = -90
    dLon_min = 180
    dLon_max = -180          


    n_colors = pLayer.GetFeatureCount()
    
    colours = cm.rainbow(np.linspace(0, 1, n_colors))
    for pFeature in pLayer:
        pGeometry_in = pFeature.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()
        if iFlag_label ==1:
            iSegment = pFeature.GetField("iseg")
        if sGeometry_type =='LINESTRING':
            dummy0 = loads( pGeometry_in.ExportToWkt() )
            aCoords_gcs = dummy0.coords
            aCoords_gcs= np.array(aCoords_gcs)
            aCoords_gcs = aCoords_gcs[:,0:2]
            nvertex = len(aCoords_gcs)
            
            for i in range(nvertex):
                dLon = aCoords_gcs[i][0]
                dLat = aCoords_gcs[i][1]
                if dLon > dLon_max:
                    dLon_max = dLon
                
                if dLon < dLon_min:
                    dLon_min = dLon
                
                if dLat > dLat_max:
                    dLat_max = dLat

                if dLat < dLat_min:
                    dLat_min = dLat
                
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
            if n_colors < 10:
                line, = ax.plot(x, y, color= colours[lID],linewidth=1)
            else:
                line, = ax.plot(x, y, color= 'black',linewidth=1)
            lID = lID + 1
            #add label 
            if iFlag_label ==1:
                sText = 'Seg ' + "{:0d}".format( iSegment )
                ax.text(dLon_label,dLat_label, sText, \
                    verticalalignment='center', horizontalalignment='center',\
                    #transform=ax.transAxes, \
                    color='black', fontsize=7)
            

    pDataset = pLayer = pFeature  = None    
   
    

    if aExtent_in is None:
        marginx  = (dLon_max - dLon_min) / 20
        marginy  = (dLat_max - dLat_min) / 20
        aExtent = [dLon_min - marginx , dLon_max + marginx , dLat_min - marginy , dLat_max + marginy]
    else:
        aExtent = aExtent_in  
 
    
    ax.set_extent(aExtent)       

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.3, linestyle='--')
    gl.xlocator = mticker.MaxNLocator(5)
    gl.ylocator = mticker.MaxNLocator(5)
    gl.xlabel_style = {'size': 8, 'color': 'k', 'rotation':0, 'ha':'right'}
    gl.ylabel_style = {'size': 8, 'color': 'k', 'rotation':90,'weight': 'normal'}
    if iFlag_title is None:
        ax.set_title( sTitle )
    else:
        if iFlag_title==1:
            ax.set_title( sTitle )
        else:
            pass
         
    
    sFilename_out = os.path.join(self.sWorkspace_output_basin, sFilename_in)
    plt.savefig(sFilename_out, bbox_inches='tight')
    plt.close(fig)
    #plt.clf()

    return

def _plot_area_of_difference(self, iCase_index, iMesh_type, sMesh_type, sFilename_in, aExtent_in = None, pProjection_map_in = None):
    #request = cimgt.OSM()
    sFilename_json = self.sFilename_area_of_difference
    sFilename_json = os.path.join(self.sWorkspace_output_basin, sFilename_json)
    
    fig = plt.figure( dpi=300)
    fig.set_figwidth( 4 )
    fig.set_figheight( 4 )
    ax = fig.add_axes([0.1, 0.15, 0.75, 0.8] , projection=pProjection_map ) #request.crs
    
    pDriver = ogr.GetDriverByName('GeoJSON')
    pDataset = pDriver.Open(sFilename_json, gdal.GA_ReadOnly)
    pLayer = pDataset.GetLayer(0)

    pSrs = osr.SpatialReference()  
    pSrs.ImportFromEPSG(4326)    # WGS84 lat/lon

    lID = 0
    dLat_min = 90
    dLat_max = -90
    dLon_min = 180
    dLon_max = -180          
    
    #ax.add_image(request, 6)    # 5 = zoom level
    n_colors = pLayer.GetFeatureCount()
    
    colours = cm.rainbow(np.linspace(0, 1, n_colors))
    lID=0
    for pFeature in pLayer:
        pGeometry_in = pFeature.GetGeometryRef()
        sGeometry_type = pGeometry_in.GetGeometryName()
        if sGeometry_type =='POLYGON':
            dummy0 = loads( pGeometry_in.ExportToWkt() )
            aCoords_gcs = dummy0.exterior.coords
            aCoords_gcs= np.array(aCoords_gcs)
            nvertex = len(aCoords_gcs)
            for i in range(nvertex):
                dLon = aCoords_gcs[i][0]
                dLat = aCoords_gcs[i][1]
                if dLon > dLon_max:
                    dLon_max = dLon
                
                if dLon < dLon_min:
                    dLon_min = dLon
                
                if dLat > dLat_max:
                    dLat_max = dLat

                if dLat < dLat_min:
                    dLat_min = dLat

            polygon = mpatches.Polygon(aCoords_gcs[:,0:2], closed=True,  linewidth=0.25, \
                alpha=0.8, edgecolor = 'red',facecolor='red', \
                    transform=ccrs.PlateCarree() )
            ax.add_patch(polygon)   
            lID = lID + 1
            

    pDataset = pLayer = pFeature  = None    
    
    
    
    if aExtent_in is None:
        marginx  = (dLon_max - dLon_min) / 20
        marginy  = (dLat_max - dLat_min) / 20
        aExtent = [dLon_min - marginx , dLon_max + marginx , dLat_min - marginy , dLat_max + marginy]
    else:
        aExtent = aExtent_in
        
    ax.set_extent( aExtent )      
          
    sTitle = 'Area of difference'
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=1, color='gray', alpha=0.3, linestyle='--')
    gl.xlabel_style = {'size': 8, 'color': 'k', 'rotation':0, 'ha':'right'}
    gl.ylabel_style = {'size': 8, 'color': 'k', 'rotation':90,'weight': 'normal'}
    ax.set_title( sTitle)  

    iCase = iCase_index
    if sMesh_type == 'hexagon':
        iCase = 10 - iCase
    if sMesh_type == 'square':
        iCase = 10 - iCase
    if sMesh_type == 'latlon':
        iCase = 10 - iCase 
    sText = 'Case: ' + "{:0d}".format( int(iCase) ) 
    ax.text(0.05, 0.95, sText, \
    verticalalignment='top', horizontalalignment='left',\
            transform=ax.transAxes, \
            color='black', fontsize=8)

    sText = 'Mesh type: ' + sMesh_type.title()
    ax.text(0.05, 0.90, sText, \
    verticalalignment='top', horizontalalignment='left',\
            transform=ax.transAxes, \
            color='black', fontsize=8)
    #we need to read the dArea_of_difference

    sFilename_basin_info = os.path.join(self.sWorkspace_output_basin, self.sFilename_basin_info)
    with open(sFilename_basin_info) as json_file:
        aBasin_info = json.load(json_file) 
        self.dArea_of_difference = float(aBasin_info['dArea_of_difference'])

    

    sText = 'Total area: ' + "{:4.1f}".format( int(self.dArea_of_difference/1.0E6) ) + ' km^2'
    ax.text(0.05, 0.85, sText, \
    verticalalignment='top', horizontalalignment='left',\
            transform=ax.transAxes, \
            color='blue', fontsize=8)
    
    sFilename_out = os.path.join(self.sWorkspace_output_basin, sFilename_in)
    plt.savefig(sFilename_out, bbox_inches='tight')
    #plt.show()
    return