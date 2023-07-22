import os
from pathlib import Path
#dependency packages
from pyflowline.external.pyearth.visual.map.map_vector_polygon_data import map_vector_polygon_data
from pyflowline.external.pyearth.visual.map.map_vector_polyline_data import map_vector_polyline_data
from pyflowline.external.pyearth.visual.map.map_multiple_vector_data import map_multiple_vector_data


def basin_plot(self,
               iFlag_type_in,
               iMesh_type,     
               sMesh_type,
               sFilename_output_in=None,
               sFilename_mesh_in = None,
               iFlag_title_in=None,
               sVariable_in=None,
               aExtent_in = None,
               pProjection_map_in = None):
    """_summary_

    Args:
        iMesh_type (_type_): _description_
        sMesh_type (_type_): _description_
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
    
    sField_thickness = ''
    iFlag_thickness = 0
    if iFlag_type_in ==1:
        #point based
        pass
    else:
        if iFlag_type_in == 2: #polyline based
            self._plot_polyline_variable(sVariable_in,
                                               sFilename_output_in=sFilename_output_in,
                                               aExtent_in = aExtent_in,
                                               pProjection_map_in = pProjection_map_in)
            pass
        else:
            if iFlag_type_in == 3:#polygon based

                self._plot_polygon_variable( sVariable_in,                 
                                               sFilename_output_in=sFilename_output_in,
                                               aExtent_in = aExtent_in,
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
                            sFilename1 = ''
                            aFiletype_in = [3, 2, 2]
                            aFilename_in = [sFilename_mesh, sFilename0, sFilename1]
                            map_multiple_vector_data(aFiletype_in,
                                             aFilename_in,
                                             sFilename_output_in=sFilename_output_in,
                                             sTitle_in= 'Mesh with flowline and observation',
                                             aFlag_color_in = [0, 0, 0],
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
                             dData_min_in = None,
                             dData_max_in = None,
                             sFilename_output_in=None,
                             aExtent_in = None,
                             pProjection_map_in = None):
    


    iFlag_label = 0
    iFlag_thickness = 0
    if sVariable_in is not None:
        if sVariable_in == 'flowline_raw':
            sFilename_json = self.sFilename_flowline_raw
            sTitle = 'Original flowline'
            iFlag_color = 0
            iFlag_thickness = 0
        else:
            if sVariable_in == 'flowline_filter':
                sFilename_json = self.sFilename_flowline_filter
                sTitle = 'Filtered flowline'
                iFlag_color = 0
                iFlag_thickness = 0
            else:
                if sVariable_in == 'flowline_simplified':
                    sFilename_json = self.sFilename_flowline_simplified                    
                    sTitle = 'Simplified flowline'
                    iFlag_color = 1
                    iFlag_thickness = 0
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
                                self._plot_area_of_difference( self.iMesh_type,
                                                         self.sMesh_type,
                                                         sFilename_output_in,
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
               
    
    map_vector_polyline_data(sFilename_json,
                             sFilename_output_in= sFilename_output_in,
                             iFlag_title_in=iFlag_title_in,
                             iFlag_thickness_in= iFlag_thickness  ,
                             sTitle_in=sTitle,
                             iFlag_color_in= iFlag_color,
                             iFlag_label_in=iFlag_label,
                             sField_thickness_in = sField_thickness,
                             aExtent_in = aExtent_in,
                             pProjection_map_in = pProjection_map_in)
    
    

def _plot_polygon_variable(self,
                             sVariable_in,
                             iFigwidth_in=None,
                             iFigheight_in=None,
                             dData_min_in = None,
                             dData_max_in = None,
                             sFilename_output_in=None,
                             aExtent_in = None,
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
            dData_min = dData_min_in
            dData_max = dData_max_in
            sFilename = self.sFilename_elevation
            sFilename = self.sFilename_variable_polygon
        else:
            if sVariable_in == 'drainage_area':
                sVariable='drainage_area'
                sTitle = 'Drainage area'
                sUnit = r'Unit: $m^{2}$'
                dData_min = dData_min_in
                dData_max = dData_max_in
                sFilename = self.sFilename_variable_polygon
            else:
                if sVariable_in == 'travel_distance':
                    sVariable='travel_distance'
                    sTitle = 'Distance to outlet'
                    sUnit = r'Unit: m'
                    dData_min = 0.0
                    dData_max = dData_max_in
                    sFilename = self.sFilename_variable_polygon
                else:
                    sVariable='slope'
                    sTitle = 'Surface slope'
                    sUnit = 'Unit: percent'
                    dData_min = 0.0
                    dData_max = dData_max_in
                    sFilename = self.sFilename_variable_polygon
        
    else:
        if sVariable_in == 'elevation':
            sVariable='elevation'
            sTitle = 'Surface elevation'
            sUnit = r'Unit: m'
            dData_min = dData_min_in
            dData_max = dData_max_in
            sFilename = self.sFilename_variable_polygon
        else:
            if sVariable_in == 'drainage_area':
                sVariable='drainage_area'
                sTitle = 'Drainage area'
                sUnit = r'Unit: $m^{2}$'
                dData_min = dData_min_in
                dData_max = dData_max_in
                sFilename = self.sFilename_variable_polygon

            else:
                if sVariable_in == 'travel_distance':
                    sVariable='travel_distance'
                    sTitle = 'Travel distance'
                    sUnit = r'Unit: m'
                    dData_min = 0.0
                    dData_max = dData_max_in
                    iFlag_subbasin = 1
                    sFilename = self.sFilename_variable_polygon
                else:
                    sVariable='slope'
                    sTitle = 'Surface slope'
                    sUnit = r'Unit: percent'
                    dData_min = dData_min_in
                    dData_max = dData_max_in
                    sFilename = self.sFilename_variable_polygon
        pass
    
    map_vector_polygon_data(sFilename,
                             iFlag_colorbar_in = 1,
                             sFilename_output_in=sFilename_output_in,
                             sTitle_in= sTitle,
                             aExtent_in = aExtent_in,
                             pProjection_map_in = pProjection_map_in)


    return


#this is a reserved function
def _plot_area_of_difference(self, iMesh_type, sMesh_type, sFilename_figure_in, aExtent_in = None, pProjection_map_in = None):

    sFilename_json = self.sFilename_area_of_difference    

    sFilename_in = self.sFilename_mesh
    sFilename_out = sFilename_figure_in

    map_vector_polygon_data(1, sFilename_in, sFilename_out, 'cell')


    return