import os
#dependency packages
from pyflowline.external.pyearth.visual.map.map_vector_polygon_data import map_vector_polygon_data
from pyflowline.external.pyearth.visual.map.map_vector_polyline_data import map_vector_polyline_data
from pyflowline.external.pyearth.visual.map.map_multiple_vector_data import map_multiple_vector_data

#plot function
def plot(self,
         iFlag_type_in = None,
         iFlag_title_in = None,
         sFilename_output_in=None,
         sVariable_in=None,
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
    

    if iFlag_type_in == 1: #point based, such as dam
        #currently, this feature is not supported
        pass
    else:
        if iFlag_type_in == 2:
            #polyline based, only flowline
            for pBasin in self.aBasin:
                pBasin.basin_plot(self.iCase_index,
                                  self.iMesh_type,
                                  self.sMesh_type,
                                  iFlag_title_in= iFlag_title_in,
                                  sFilename_output_in=sFilename_output_in,
                                  sVariable_in= sVariable_in,
                                  aExtent_in=aExtent_in,
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

def _plot_mesh(self, sFilename_output_in=None, aExtent_in=None, pProjection_map_in = None):

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
                             iFlag_title=None,
                             aExtent_in=None,
                             pProjection_map_in = None):

    aFiletype_in = list()
    aFilename_in = list()
    aFilename_in.append(self.sFilename_mesh)
    aFiletype_in.append(1)

    for pBasin in self.aBasin:
        aFiletype_in.append(2)
        dummy = pBasin.sFilename_flowline_conceptual
        sFilename_json = os.path.join(pBasin.sWorkspace_output_basin, dummy)
        aFilename_in.append(sFilename_json)

    map_multiple_vector_data(aFiletype_in,
                             aFilename_in,
                             sFilename_output_in=sFilename_output_in,
                             sTitle_in= 'Mesh with flowline',
                             aFlag_color_in=[0, 1],
                             aExtent_in = aExtent_in,
                             pProjection_map_in = pProjection_map_in)
    return

def basin_plot(self,
               iCase_index,
               iMesh_type,
               sMesh_type,
               sFilename_output_in=None,
               iFlag_title_in=None,
               sVariable_in=None,
               aExtent_in = None,
               pProjection_map_in = None):
    """_summary_

    Args:
        iCase_index (_type_): _description_
        iMesh_type (_type_): _description_
        sMesh_type (_type_): _description_
        sFilename_output_in (_type_, optional): _description_. Defaults to None.
        iFlag_title_in (_type_, optional): _description_. Defaults to None.
        sVariable_in (_type_, optional): _description_. Defaults to None.
        aExtent_in (_type_, optional): _description_. Defaults to None.
        pProjection_map_in (_type_, optional): _description_. Defaults to None.
    """


    iFlag_label = 0
    sWorkspace_output_basin = self.sWorkspace_output_basin
    if sVariable_in is not None:
        if sVariable_in == 'flowline_raw':
            sFilename_json = self.sFilename_flowline_raw
            sTitle = 'Original flowline'
            iFlag_color = 0
        else:
            if sVariable_in == 'flowline_filter':
                sFilename_json = self.sFilename_flowline_filter
                sTitle = 'Filtered flowline'
                iFlag_color = 0
            else:
                if sVariable_in == 'flowline_simplified':
                    sFilename_json = self.sFilename_flowline_simplified                    
                    sTitle = 'Simplified flowline'
                    iFlag_color = 1
                    if aExtent_in is None:
                        iFlag_label = 1
                    else:
                        iFlag_label=0
                else:
                    if sVariable_in == 'flowline_conceptual':
                        sFilename_json = self.sFilename_flowline_conceptual                        
                        sTitle = 'Conceptual flowline'
                        iFlag_color = 1
                        if aExtent_in is None:
                            iFlag_label = 1
                        else:
                            iFlag_label=0
                    else:
                        if sVariable_in == 'aof':
                            sFilename_json = self.sFilename_area_of_difference                            
                            sTitle = 'Conceptual flowline'
                            iFlag_label = 1
                            self.plot_area_of_difference(iCase_index,
                                                         iMesh_type,
                                                         sMesh_type,
                                                         sFilename_output_in,
                                                         aExtent_in = aExtent_in)
                            return
                        else:
                            pass
                pass
    else:
        #default
        sFilename_json = self.sFilename_flowline_conceptual


    map_vector_polyline_data(sFilename_json,
                             sFilename_output_in,
                             iFlag_title_in=iFlag_title_in,
                             iFlag_thickness_in=0  ,
                             sTitle_in=sTitle,
                             iFlag_color_in= iFlag_color,
                             iFlag_label_in=iFlag_label,
                             aExtent_in = aExtent_in,
                             pProjection_map_in = pProjection_map_in)

    return

#this function is used to plot the study area
def _plot_study_area(self, sFilename_boundary_in = None, sFilename_slope_in = None, sFilename_nhd_in = None):
    return

#this is a reserved function
def _compare_with_raster_dem_method(self, sFilename_dem_flowline, sFilename_in, aExtent_in=None, pProjection_map_in = None):

    return

#this is a reserved function
def _plot_area_of_difference(self, iCase_index, iMesh_type, sMesh_type, sFilename_figure_in, aExtent_in = None, pProjection_map_in = None):

    sFilename_json = self.sFilename_area_of_difference    

    sFilename_in = self.sFilename_mesh
    sFilename_out = sFilename_figure_in

    map_vector_polygon_data(1, sFilename_in, sFilename_out, 'cell')


    return
