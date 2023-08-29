import os
from pathlib import Path
#dependency packages
from pyflowline.external.pyearth.visual.map.map_vector_polygon_data import map_vector_polygon_data
from pyflowline.external.pyearth.visual.map.map_vector_polyline_data import map_vector_polyline_data
from pyflowline.external.pyearth.visual.map.map_multiple_vector_data import map_multiple_vector_data

def replace_last_occurrence(sFilename_path_in, sSubstring_in, sSubstring_out):
    last_occurrence_index = sFilename_path_in.rfind(sSubstring_in)
    if last_occurrence_index == -1:
        # Substring not found, return the original string
        return sFilename_path_in
    else:
        return sFilename_path_in[:last_occurrence_index] + sSubstring_out + sFilename_path_in[last_occurrence_index+len(sSubstring_in):]



def basin_plot(self,
               iFlag_type_in,               
               sMesh_type,
               sFilename_output_in=None,
               sFilename_mesh_in = None,
               iFont_size_in = None,
               iFlag_title_in=None,
               iFlag_colorbar_in=None,
               iFlag_scientific_notation_colorbar_in=None,
               dData_min_in = None,
               dData_max_in = None,
               sVariable_in=None,
               aExtent_in = None,
                aLegend_in = None,
               pProjection_map_in = None):
    """_summary_

    Args:
        sMesh_type (_type_): For labeling purpose
        sFilename_output_in (_type_, optional): _description_. Defaults to None.
        iFlag_title_in (_type_, optional): _description_. Defaults to None.
        sVariable_in (_type_, optional): _description_. Defaults to None.
        aExtent_in (_type_, optional): _description_. Defaults to None.
        pProjection_map_in (_type_, optional): _description_. Defaults to None.
    """

    if sFilename_mesh_in is None:
        sFilename_mesh = self.sFilename_mesh
    else:
        sFilename_mesh = sFilename_mesh_in
    

    if iFlag_type_in ==1:
        #point based
        pass
    else:
        if iFlag_type_in == 2: #polyline based
            self._plot_polyline_variable(sVariable_in,
                                         iFlag_title_in=iFlag_title_in,
                                           iFont_size_in=iFont_size_in,
                                               sFilename_output_in=sFilename_output_in,
                                               aExtent_in = aExtent_in,
                                               aLegend_in = aLegend_in,
                                               pProjection_map_in = pProjection_map_in)
            pass
        else:
            if iFlag_type_in == 3:#polygon based

                self._plot_polygon_variable( sVariable_in,                 
                                            iFlag_title_in= iFlag_title_in,
                                            iFont_size_in=iFont_size_in,
                                            iFlag_colorbar_in=iFlag_colorbar_in,
                                            iFlag_scientific_notation_colorbar_in=iFlag_scientific_notation_colorbar_in,
                                            dData_min_in = dData_min_in,
                                            dData_max_in = dData_max_in,
                                               sFilename_output_in=sFilename_output_in,
                                               aExtent_in = aExtent_in,
                                                aLegend_in = aLegend_in,
                                               pProjection_map_in = pProjection_map_in)

                pass
            else:
                if iFlag_type_in == 4:#mixed
                    if sVariable_in == "flow_direction_with_mesh":
                        sFilename = self.sFilename_flow_direction #this can be either domain wide or subbasin level
                        aFiletype_in = [3, 2]
                        aFilename_in = [sFilename_mesh, sFilename]
                        map_multiple_vector_data(aFiletype_in,
                                             aFilename_in,
                                             sFilename_output_in=sFilename_output_in,
                                             sTitle_in= 'Mesh with flowline',
                                             aFlag_color_in=[0, 0],
                                             aFlag_fill_in = [0,0])
                    else:
                        if sVariable_in == "flow_direction_with_observation":
                            sFilename0 = self.sFilename_flow_direction #this can be either domain wide or subbasin level
                            #should use the pyflowline simplified flowline
                            sFilename_dummy = self.sFilename_flowline_simplified
                            #now replace the folder string 
                            sFilename1 = replace_last_occurrence(sFilename_dummy, 'hexwatershed', 'pyflowline')
                            sFilename_dummy = self.sFilename_flowline_conceptual
                            #now replace the folder string 
                            sFilename2 = replace_last_occurrence(sFilename_dummy, 'hexwatershed', 'pyflowline')

                            aFiletype_in = [2, 2, 2]
                            aFilename_in = [sFilename0, sFilename1, sFilename2]
                            map_multiple_vector_data(aFiletype_in,
                                             aFilename_in,                           
                                             iFont_size_in=iFont_size_in,                 
                                             sFilename_output_in=sFilename_output_in,
                                             sTitle_in= 'Flow direction with observation',
                                             aFlag_thickness_in=  [1, 0, 0],
                                             aVariable_in= ['drainage_area', '', 'stream_segment'],
                                             aLegend_in = aLegend_in,
                                             aFlag_color_in = [0, 0, 1],
                                             aFlag_fill_in  = [0, 0, 0])
                        else:
                            print('Unsupported variable: ', sVariable_in, ' in basin_plot.')
                            return                    
                    
                    pass
                else:
                    #unsupported
                    pass   

    return

def _plot_polyline_variable(self,
                             sVariable_in,
                             iFlag_title_in=None,
                             iFigwidth_in=None,
                             iFigheight_in=None,
                             iFont_size_in=None,
                             dData_min_in = None,
                             dData_max_in = None,
                             sFilename_output_in=None,
                             aExtent_in = None,
                             aLegend_in = None,
                             pProjection_map_in = None):
    
    


    iFlag_label = 0
    iFlag_thickness = 0
    if sVariable_in is not None:
        if sVariable_in == 'flowline_raw':
            sFilename_json = self.sFilename_flowline_raw
            sTitle = 'Original flowline'
            iFlag_color = 0
            iFlag_thickness = 0
            sField_thickness = None
        else:
            if sVariable_in == 'flowline_filter':
                sFilename_json = self.sFilename_flowline_filter
                sTitle = 'Filtered flowline'
                iFlag_color = 0                
                iFlag_thickness = 0
                sField_thickness = None
            else:
                if sVariable_in == 'flowline_simplified':
                    sFilename_json = self.sFilename_flowline_simplified                    
                    sTitle = 'Simplified flowline'
                    iFlag_color = 1
                    iFlag_thickness = 0
                    sField_thickness = None

                    if aExtent_in is None:
                        iFlag_label = 1
                    else:
                        iFlag_label=0
                else:
                    if sVariable_in == 'flowline_conceptual':
                        sFilename_json = self.sFilename_flowline_conceptual                        
                        sTitle = 'Conceptual flowline'
                        iFlag_color = 1
                        iFlag_thickness = 0
                        sField_thickness = None

                        if aExtent_in is None:
                            iFlag_label = 1
                        else:
                            iFlag_label=0
                    else:
                        if sVariable_in == 'flow_direction':
                            sFilename_json = self.sFilename_flow_direction
                            iFlag_label= 0
                            iFlag_color = 0 
                            iFlag_thickness = 1                                    
                            sField_thickness = 'drainage_area'
                            sTitle = 'Flow direction'                                               
                        else:
                            if sVariable_in ==  'aof':
                                sFilename_json = self.sFilename_area_of_difference                            
                                sTitle = 'Conceptual flowline'
                                iFlag_label = 1
                                self._plot_area_of_difference( sFilename_output_in,
                                                         aExtent_in = aExtent_in)
                                return
                            else:
                                print('Unsupported variable: ', sVariable_in, ' in basin_plot.')
                            pass
                pass
    else:
        #default
        print('A variable is needed.')
        return
    
    if iFlag_title_in is not None:        
        if iFlag_title_in == 0:
            sTitle=''
        else:
            pass
    else:        
        sTitle=''
        pass
               
    
    map_vector_polyline_data(sFilename_json,
                             sFilename_output_in= sFilename_output_in,                             
                             iFlag_thickness_in= iFlag_thickness  ,
                             sTitle_in=sTitle,
                             iFlag_color_in= iFlag_color,
                             iFlag_label_in=iFlag_label,
                              iFont_size_in=iFont_size_in,
                             sField_thickness_in = sField_thickness,
                             aExtent_in = aExtent_in,
                             aLegend_in = aLegend_in,
                             pProjection_map_in = pProjection_map_in)
    
def _plot_polygon_variable(self,
                             sVariable_in,
                             iFigwidth_in=None,
                             iFigheight_in=None,
                             iFlag_colorbar_in=None,
                            iFlag_title_in=None,
                            iFont_size_in=None,
                            iFlag_scientific_notation_colorbar_in = None,
                             dData_min_in = None,
                             dData_max_in = None,
                             sFilename_output_in=None,
                              aExtent_in = None,
                             aLegend_in = None,
                             pProjection_map_in = None):
    """_summary_

    Args:
        sVariable_in (_type_): _description_
        sFilename_output_in (_type_, optional): _description_. Defaults to None.
        iFigwidth_in (_type_, optional): _description_. Defaults to None.
        iFigheight_in (_type_, optional): _description_. Defaults to None.
        aExtent_in (_type_, optional): _description_. Defaults to None.
        pProjection_map_in (_type_, optional): _description_. Defaults to None.
        dData_min_in (_type_, optional): _description_. Defaults to None.
        dData_max_in (_type_, optional): _description_. Defaults to None.
    """
    
    sMesh_type = self.sMesh_type
   
    if sMesh_type == 'mpas':
        if sVariable_in == 'elevation':
            sVariable='elevation' #Elevation_profile'
            sTitle = 'Surface elevation'
            sUnit = 'Unit: m'
            sColormap ='terrain'
            dData_min = dData_min_in
            dData_max = dData_max_in
            sFilename = self.sFilename_elevation
            sFilename = self.sFilename_variable_polygon
        else:
            if sVariable_in == 'drainage_area':
                sVariable='drainage_area'
                sTitle = 'Drainage area'
                sUnit = r'Units: $m^{2}$'
                dData_min = dData_min_in
                dData_max = dData_max_in
                sColormap ='Spectral_r'
                sFilename = self.sFilename_variable_polygon
            else:
                if sVariable_in == 'travel_distance':
                    sVariable='travel_distance'
                    sTitle = 'Distance to outlet'
                    sUnit = r'Unit: m'
                    dData_min = 0.0
                    dData_max = dData_max_in
                    sColormap ='Spectral_r'
                    sFilename = self.sFilename_variable_polygon
                else:
                    sVariable='slope'
                    sTitle = 'Surface slope'
                    sUnit = 'Unit: percent'
                    sColormap='Spectral_r'
                    dData_min = 0.0
                    dData_max = dData_max_in
                    sFilename = self.sFilename_variable_polygon
        
    else:
        if sVariable_in == 'area':
            sVariable='area'
            sTitle = 'Area'
            sUnit = r'Units: $m^{2}$'
            sColormap ='terrain'
            dData_min = dData_min_in
            dData_max = dData_max_in
            sFilename = self.sFilename_variable_polygon
            pass
        else:
            if sVariable_in == 'elevation':
                sVariable='elevation'
                sTitle = 'Surface elevation'
                sUnit = r'Unit: m'
                sColormap ='terrain'
                dData_min = dData_min_in
                dData_max = dData_max_in
                sFilename = self.sFilename_variable_polygon
            else:
                if sVariable_in == 'drainage_area':
                    sVariable='drainage_area'
                    sTitle = 'Drainage area'
                    sUnit = r'Units: $m^{2}$'
                    dData_min = dData_min_in
                    dData_max = dData_max_in
                    sColormap ='Spectral_r'
                    sFilename = self.sFilename_variable_polygon

                else:
                    if sVariable_in == 'travel_distance':
                        sVariable='travel_distance'
                        sTitle = 'Travel distance'
                        sUnit = r'Unit: m'
                        dData_min = 0.0
                        dData_max = dData_max_in
                        sColormap ='Spectral_r'
                        iFlag_subbasin = 1
                        sFilename = self.sFilename_variable_polygon
                    else:
                        sVariable='slope'
                        sTitle = 'Surface slope'
                        sUnit = r'Unit: percent'
                        sColormap='Spectral_r'
                        dData_min = dData_min_in
                        dData_max = dData_max_in
                        sFilename = self.sFilename_variable_polygon
            pass
    
    if iFlag_title_in is not None:        
        if iFlag_title_in == 0:
            sTitle=''
        else:
            pass
    else:        
        sTitle=''
        pass

    map_vector_polygon_data(sFilename,
                            iFlag_color_in = 1,
                             iFlag_colorbar_in = iFlag_colorbar_in,
                             iFont_size_in = iFont_size_in,
                             iFlag_scientific_notation_colorbar_in = iFlag_scientific_notation_colorbar_in,
                             dData_max_in = dData_max,
                             dData_min_in = dData_min,
                             sFilename_output_in=sFilename_output_in,
                             sVariable_in= sVariable,
                             sTitle_in= sTitle,
                             sUnit_in= sUnit,
                             sColormap_in = sColormap,
                             aExtent_in = aExtent_in,
                             aLegend_in = aLegend_in,
                             pProjection_map_in = pProjection_map_in)


    return

def _plot_mesh_with_flowline(self,
                             sFilename_output_in=None,
                             iFlag_title_in=None,
                             aExtent_in=None,
                             pProjection_map_in = None):

    aFiletype_in = list()
    aFilename_in = list()
    aFlag_color = list()
    aFilename_in.append(self.sFilename_mesh)
    aFiletype_in.append(3)
    aFlag_color.append(0)

    
    dummy = self.sFilename_flowline_conceptual
    sFilename_json = os.path.join(self.sWorkspace_output_basin, dummy)
    aFilename_in.append(sFilename_json)
    aFlag_color.append(1)

    map_multiple_vector_data(aFiletype_in,
                             aFilename_in,
                             sFilename_output_in=sFilename_output_in,
                             sTitle_in= 'Mesh with flowline',
                             aFlag_color_in=aFlag_color,
                             aExtent_in = aExtent_in,
                             pProjection_map_in = pProjection_map_in)
    return

def _plot_mesh_with_flow_direction(self,
                                   iFigwidth_in=None,
                                   iFigheight_in=None,
                                   sMesh_type_in = None,
                                   sFilename_mesh_in = None,
                                   sFilename_output_in = None,
                                   aExtent_in = None,
                                   pProjection_map_in = None):
    if sMesh_type_in is None:
        sMesh_type = self.sMesh_type
    else:
        sMesh_type = sMesh_type_in

    if sFilename_mesh_in is None:
        sFilename_mesh = self.sFilename_mesh
    else:
        sFilename_mesh = sFilename_mesh_in

    sFilename = self.sFilename_flow_direction #this can be either domain wide or subbasin level

    aFiletype_in = [3, 2]

    aFilename_in = [sFilename_mesh, sFilename]
    map_multiple_vector_data(aFiletype_in,
                             aFilename_in,
                             sFilename_output_in=sFilename_output_in,
                             sTitle_in= 'Mesh with flowline',
                             aFlag_color_in=[0, 1])
    return

#this is a reserved function
def _plot_area_of_difference(self, sFilename_figure_in, aExtent_in = None, pProjection_map_in = None):

    sFilename_json = self.sFilename_area_of_difference    

    sFilename_in = self.sFilename_mesh
    sFilename_out = sFilename_figure_in

    map_vector_polygon_data(1, sFilename_in, sFilename_out, 'cell')


    return
