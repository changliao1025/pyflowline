import os

#dependency packages
from pyearth.visual.map.vector.map_vector_polygon_data import map_vector_polygon_data
from pyearth.visual.map.vector.map_vector_polyline_data import map_vector_polyline_data
from pyearth.visual.map.vector.map_multiple_vector_data import map_multiple_vector_data

#plot function

def plot(self,
         iFlag_type_in = None,
         iFlag_title_in = None,
         sFilename_output_in = None,
         sVariable_in = None,
         aExtent_in = None,
         pProjection_map_in = None):
    """_summary_

    Args:
        iFlag_type_in (_type_, optional): _description_. Defaults to None.
        iFlag_title_in (_type_, optional): _description_. Defaults to None.
        sFilename_output_in (_type_, optional): _description_. Defaults to None.
        sVariable_in (_type_, optional): _description_. Defaults to None.
        aExtent_in (_type_, optional): _description_. Defaults to None.
        pProjection_map_in (_type_, optional): _description_. Defaults to None.
    """

    if iFlag_type_in is None:
        iFlag_type_in = 2 #polyline based, only flowline

    if iFlag_title_in is None:
        iFlag_title_in = 1

    if sVariable_in is None:
        sVariable_in = 'flowline_conceptual'
    else:
        if sVariable_in == 'mesh': 
            iFlag_type_in = 3
        else:
            if sVariable_in == 'overlap':
                iFlag_type_in = 4
    

    if iFlag_type_in == 1: #point based, such as dam
        #currently, this feature is not supported
        pass
    else:
        if iFlag_type_in == 2:
            #polyline based, only flowline
            aLegend = list()
            sText = 'Case: ' + "{:0d}".format( int(self.iCase_index) ) 
            aLegend.append(sText)
            sText = 'Mesh type: ' + self.sMesh_type.title()
            aLegend.append(sText)
            sResolution =  'Resolution: ' + "{:0d}".format( int(self.dResolution_meter) ) + 'm'
            aLegend.append(sResolution)            
            for pBasin in self.aBasin:                
                pBasin.basin_plot(iFlag_type_in,
                                  self.sMesh_type,
                                  iFlag_title_in= iFlag_title_in,
                                  sFilename_output_in=sFilename_output_in,
                                  sVariable_in= sVariable_in,
                                  aExtent_in = aExtent_in,
                                  aLegend_in = aLegend,
                                  pProjection_map_in = pProjection_map_in)
        else:
            if iFlag_type_in == 3: #polygon based
                if sVariable_in == 'mesh':
                    self._plot_mesh(sFilename_output_in=sFilename_output_in,
                                    aExtent_in = aExtent_in,
                                    pProjection_map_in = pProjection_map_in)

            else:
                if iFlag_type_in == 4: #mixed based mesh + flowline
                    if sVariable_in == 'overlap':
                        self._plot_mesh_with_flowline(sFilename_output_in=sFilename_output_in,
                                                      iFlag_title_in= iFlag_title_in,
                                                      aExtent_in=aExtent_in,
                                                      pProjection_map_in = pProjection_map_in)



                    pass
                else:
                    pass

        return

def _plot_mesh(self, 
               sFilename_output_in=None, 
               aExtent_in=None, 
               pProjection_map_in = None):

    sFilename_in = self.sFilename_mesh 
    sMesh_type = self.sMesh_type

    map_vector_polygon_data(sFilename_in,
                            sFilename_output_in = sFilename_output_in,
                            sTitle_in = sMesh_type,
                            aExtent_in = aExtent_in,
                            pProjection_map_in = pProjection_map_in)

    return

#plot both polygon and polyline
def _plot_mesh_with_flowline(self,
                             sFilename_output_in=None,
                             iFlag_title_in=None,
                             aExtent_in=None,
                             pProjection_map_in = None):

    aFiletype_in = list()
    aFilename_in = list()
    aFlag_color = list()
    aVariable_in = list()
    aFilename_in.append(self.sFilename_mesh)
    aFiletype_in.append(3)
    aFlag_color.append(0)
    
    aVariable_in.append(None)

    for pBasin in self.aBasin:
        aFiletype_in.append(2)
        dummy = pBasin.sFilename_flowline_conceptual
        sFilename_json = os.path.join(pBasin.sWorkspace_output_basin, dummy)
        aFilename_in.append(sFilename_json)
        aFlag_color.append(1)
        aVariable_in.append('stream_segment')

    map_multiple_vector_data(aFiletype_in,
                             aFilename_in,
                             sFilename_output_in=sFilename_output_in,
                             sTitle_in= 'Mesh with flowline',
                             aFlag_color_in = aFlag_color,
                             aExtent_in = aExtent_in,
                             pProjection_map_in = pProjection_map_in)
    return

#this function is used to plot the study area
def _plot_study_area(self, sFilename_boundary_in = None, sFilename_slope_in = None, sFilename_nhd_in = None):
    return

#this is a reserved function
def _compare_with_raster_dem_method(self, sFilename_dem_flowline, sFilename_in, aExtent_in=None, pProjection_map_in = None):

    return

