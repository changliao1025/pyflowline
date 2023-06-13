import os
import json
#dependency packages
import numpy as np

from pyflowline.external.pyearth.toolbox.reader.text_reader_string import text_reader_string

from pyflowline.external.pyearth.visual.map.map_vector_polygon_data import map_vector_polygon_data
from pyflowline.external.pyearth.visual.map.map_vector_polyline_data import map_vector_polyline_data
from pyflowline.external.pyearth.visual.map.map_multiple_vector_data import map_multiple_vector_data

def plot(self, sFilename_output_in=None,iFlag_title = None, sVariable_in=None, aExtent_in = None):
    if sVariable_in == 'mesh':
        self._plot_mesh(sFilename_output_in=sFilename_output_in, aExtent_in = aExtent_in)
    else:            
        if sVariable_in == 'overlap':
            self._plot_mesh_with_flowline( sFilename_output_in=sFilename_output_in, iFlag_title= iFlag_title, aExtent_in=aExtent_in)
        else:            
            for pBasin in self.aBasin:  
                pBasin.basin_plot(self.iCase_index, 
                                  self.iMesh_type, 
                                  self.sMesh_type, 
                                  sFilename_output_in=sFilename_output_in,
                    iFlag_title= iFlag_title, 
                    sVariable_in= sVariable_in, 
                    aExtent_in=aExtent_in)
            pass
    
    return

#this function is used to plot the study area
def _plot_study_area(self, sFilename_boundary_in = None, sFilename_slope_in = None, sFilename_nhd_in = None):
    return

def _plot_mesh(self, sFilename_output_in=None, aExtent_in=None, pProjection_map_in = None):

    sFilename_in = self.sFilename_mesh

    map_vector_polygon_data(sFilename_in, sFilename_output_in=sFilename_output_in)
    
      
    return


#plot both polygon and polyline
def _plot_mesh_with_flowline(self, sFilename_output_in=None, iFlag_title=None, aExtent_in=None, pProjection_map_in = None):
    
    aFiletype_in = list()
    aFilename_in = list()
    aFilename_in.append(self.sFilename_mesh)
    aFiletype_in.append(1)

    for pBasin in self.aBasin: 
        aFiletype_in.append(2)
        dummy = pBasin.sFilename_flowline_conceptual
        sFilename_json = os.path.join(pBasin.sWorkspace_output_basin, dummy)
                       
        aFilename_in.append(sFilename_json)

    
  
    map_multiple_vector_data(aFiletype_in, aFilename_in, sFilename_output_in=sFilename_output_in)
    return

#this is a reserved function
def _compare_with_raster_dem_method(self, sFilename_dem_flowline, sFilename_in, aExtent_in=None, pProjection_map_in = None):
    
    return
  
def basin_plot(self, 
               iCase_index, 
               iMesh_type, 
               sMesh_type, 
               sFilename_output_in=None, 
               iFlag_title=None, 
               sVariable_in=None, 
               aExtent_in = None, 
               pProjection_map_in = None):
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
                    if aExtent_in is None:
                        iFlag_label = 1
                    else:
                        iFlag_label=0
                else:
                    if sVariable_in == 'flowline_conceptual':
                        sFilename_out = self.sFilename_flowline_conceptual
                        sFilename_json = os.path.join(sWorkspace_output_basin, sFilename_out)
                        sTitle = 'Conceptual flowline'
                        if aExtent_in is None:
                            iFlag_label = 1
                        else:
                            iFlag_label=0
                    else:
                        if sVariable_in == 'aof':
                            sFilename_out = 'area_of_difference.geojson'
                            sFilename_json = os.path.join(sWorkspace_output_basin, sFilename_out)
                            sTitle = 'Conceptual flowline'
                            iFlag_label = 1
                            self.plot_area_of_difference(iCase_index, iMesh_type, sMesh_type, sFilename_output_in, aExtent_in = aExtent_in)
                            return
                        else:
                            pass
                pass
    else:
        #default 
        sFilename_json = self.sFilename_flowline_conceptual
    
    
    map_vector_polyline_data( sFilename_json, sFilename_output_in, iFlag_thickness_in=0        )

    return



#this is a reserved function
def _plot_area_of_difference(self, iCase_index, iMesh_type, sMesh_type, sFilename_figure_in, aExtent_in = None, pProjection_map_in = None):
    
    sFilename_json = self.sFilename_area_of_difference
    sFilename_json = os.path.join(self.sWorkspace_output_basin, sFilename_json)

    sFilename_in = self.sFilename_mesh
    sFilename_out = sFilename_figure_in

    map_vector_polygon_data(1, sFilename_in, sFilename_out, 'cell')
    
    
    return