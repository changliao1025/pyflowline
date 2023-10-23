
import os, stat
import numpy as np
from pathlib import Path
import subprocess
import datetime
from shutil import copy2
from osgeo import osr, ogr, gdal
#from shapely.wkt import loads
from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_cell

from pyflowline.external.pyearth.system.define_global_variables import *
from pyflowline.formats.convert_attributes import convert_gcs_attributes_to_cell
from pyflowline.external.pyearth.gis.gdal.gdal_functions import get_geometry_coords

pDate = datetime.datetime.today()
sDate_default = "{:04d}".format(pDate.year) + "{:02d}".format(pDate.month) + "{:02d}".format(pDate.day)

#setup common resolution
aISEA3H = [4320.490,
           2539.690,
           1480.02, 
           855.419,
           494.959,
           285.6520,
           165.058, 
           95.2636, 
           55.0226, 
           31.7596, 
           18.341,
             10.5871, 
             6.11367, 
             3.52911]

aISEA4H = [3764.92,
           1913.88,
           961.978,
           481.7710,
           241.0470,
           120.56, 
           60.2893, 
           30.147, 
           15.0741, 
           7.53719, 
           3.76863]

def generate_bash_script(sWorkspace_output):
    sName  = 'dggrid.ini'
    sFilename_configuration  =  os.path.join(sWorkspace_output,  sName )
    os.chdir(sWorkspace_output)
    #detemine the system platform
    # Determine the appropriate executable name for the platform
    system = platform.system()
    if platform.system() == 'Windows':
        sFilename_executable = 'dggrid.exe'
        iFlag_unix = 0
    else:
        sFilename_executable = './dggrid'

    if system == 'Windows':
        # execute binary on Windows
        iFlag_unix = 0
    elif system == 'Linux':
        # execute binary on Linux
        iFlag_unix = 1
    elif system == 'Darwin':
        # execute binary on macOS
        iFlag_unix = 1
    else:
        # unsupported operating system
        print('Unsupported operating system: ' + system)
        print('Please reach out to the developers for assistance.')
        #generate the bash/batch script
    if iFlag_unix == 1 :
        sFilename_bash = os.path.join(str(Path(  sWorkspace_output)  ) ,  "run_dggrid.sh" )
        ofs = open(sFilename_bash, 'w')
        sLine = '#!/bin/bash\n'
        ofs.write(sLine)
        sLine = 'cd ' +   sWorkspace_output+ '\n'
        ofs.write(sLine)
        sLine = sFilename_executable + ' ' + sFilename_configuration + '\n'
        ofs.write(sLine)
        ofs.close()
        os.chmod(sFilename_bash, stat.S_IRWXU )
    else:
        sFilename_bash = os.path.join(str(Path(  sWorkspace_output)  ) ,  "run_dggrid.bat" )
        ofs = open(sFilename_bash, 'w')
        sLine = 'cd ' +   sWorkspace_output+ '\n'
        ofs.write(sLine)
        sLine = sFilename_executable + ' ' + sFilename_configuration + '\n'
        ofs.write(sLine)
        ofs.close()
        os.chmod(sFilename_bash, stat.S_IRWXU )

    return

def find_number_range(number, aArray):
    nPoint = len(aArray)
    for i in range(nPoint):
        if aArray[i] <= number <= aArray[i+1]:
            return i  # Return the index of the range where the number falls
    return -1  


def dggrid_find_index_by_resolution(sDggrid_type, dResolution):
    if sDggrid_type == 'ISEA3H':
        #unit km
        index = find_number_range(dResolution * 0.001, aISEA3H)
        iResolution_index = index+1
        pass
    else:
        if sDggrid_type == 'ISEA4H':
            index = find_number_range(dResolution * 0.001, aISEA4H)
            iResolution_index = index + 1
            pass
        pass

    return iResolution_index

def dggrid_find_resolution_by_index(sDggrid_type, iResolution_index):
    if sDggrid_type == 'ISEA3H':
       
        dResolution = aISEA3H[iResolution_index-1 ] * 1000
        pass
    else:
        if sDggrid_type == 'ISEA4H':
            dResolution = aISEA4H[iResolution_index-1 ] * 1000
            pass
        pass

    return dResolution

def convert_dggrid_mesh_to_pyflowline_mesh(sFilename_dggrid_mesh, sFilename_mesh_pyflowline):
    
    iReturn_code = 1
    if os.path.isfile(sFilename_dggrid_mesh):
        print(sFilename_dggrid_mesh)
        pass
    else:
        print('This mesh file does not exist: ', sFilename_dggrid_mesh )
        iReturn_code = 0
        return iReturn_code
    
    if os.path.isfile(sFilename_mesh_pyflowline):
        print('This mesh file already exists: ', sFilename_mesh_pyflowline )
        os.remove(sFilename_mesh_pyflowline)

    aDggrid=list()
    aDggrid_dict = dict()
    lCellIndex = 0

    pDriver_geojson = ogr.GetDriverByName('GeoJSON')    
    pDataset_mesh = pDriver_geojson.Open(sFilename_dggrid_mesh, gdal.GA_ReadOnly)
    pLayer_mesh = pDataset_mesh.GetLayer(0)
    pSpatial_reference_out = pLayer_mesh.GetSpatialRef()
    ldefn = pLayer_mesh.GetLayerDefn()

    pSpatial_reference_gcs = osr.SpatialReference()  
    pSpatial_reference_gcs.ImportFromEPSG(4326)    # WGS84 lat/lon  
    pDataset = pDriver_geojson.CreateDataSource(sFilename_mesh_pyflowline)
    pLayer = pDataset.CreateLayer('cell', pSpatial_reference_gcs, ogr.wkbPolygon)
    # Add one attribute
    pLayer.CreateField(ogr.FieldDefn('cellid', ogr.OFTInteger64)) #long type for high resolution
    pLayer.CreateField(ogr.FieldDefn('longitude', ogr.OFTReal)) #long type for high resolution
    pLayer.CreateField(ogr.FieldDefn('latitude', ogr.OFTReal)) #long type for high resolution
    pArea_field = ogr.FieldDefn('area', ogr.OFTReal)
    pArea_field.SetWidth(20)
    pArea_field.SetPrecision(2)
    pLayer.CreateField(pArea_field)
    pLayerDefn = pLayer.GetLayerDefn()
    pFeature = ogr.Feature(pLayerDefn)   
    
    #we also need to spatial reference
    for pFeature_mesh in pLayer_mesh:
        pGeometry_mesh = pFeature_mesh.GetGeometryRef()        
        #dummy0 = loads( pGeometry_mesh.ExportToWkt() )
        #aCoords_gcs = dummy0.exterior.coords
        #aCoords_gcs= np.array(aCoords_gcs)  
        aCoords_gcs = get_geometry_coords(pGeometry_mesh)   
        dLongitude_center = np.mean(aCoords_gcs[:-1,0])
        dLatitude_center = np.mean(aCoords_gcs[:-1,1])   

        lCellID = int(pFeature_mesh.GetField("name") )
        pdggrid = convert_gcs_coordinates_to_cell(5, dLongitude_center, dLatitude_center, aCoords_gcs)  
        dArea = pdggrid.calculate_cell_area()
        pdggrid.calculate_edge_length()
        pdggrid.dLength_flowline = pdggrid.dLength #Default
        pdggrid.lCellID = lCellID      

        aNeighbor = pFeature_mesh.GetField("neighbors")   
        pdggrid.nNeighbor = len(aNeighbor)
        pdggrid.aNeighbor = list()
        for i in range(len(aNeighbor)):
            pdggrid.aNeighbor.append( int(aNeighbor[i]) )       

        aDggrid.append(pdggrid)
        aDggrid_dict[lCellID] = lCellIndex
        lCellIndex = lCellIndex + 1

        nVertex = pdggrid.nVertex 
        ring = ogr.Geometry(ogr.wkbLinearRing)        
        for j in range(nVertex):
            x1 = pdggrid.aVertex[j].dLongitude_degree
            y1 = pdggrid.aVertex[j].dLatitude_degree
            ring.AddPoint(x1, y1)           
            pass
        x1 = pdggrid.aVertex[0].dLongitude_degree
        y1 = pdggrid.aVertex[0].dLatitude_degree
        ring.AddPoint(x1, y1) #double check            
        pPolygon = ogr.Geometry(ogr.wkbPolygon)
        pPolygon.AddGeometry(ring)
        pFeature.SetGeometry(pPolygon)
        pFeature.SetField("cellid", int(lCellID) )
        pFeature.SetField("longitude", dLongitude_center )
        pFeature.SetField("latitude", dLatitude_center )
        pFeature.SetField("area", dArea )
      
        pLayer.CreateFeature(pFeature)

    #rebuild neighbor list    
    aDggrid_out = list()
    
    ncell = len(aDggrid)
    aCellID  = list()
 
    for pCell in aDggrid:           
        aNeighbor = pCell.aNeighbor       
        aNeighbor_new = list()       
        for lNeighbor in aNeighbor:            
            if lNeighbor in aDggrid_dict:                
                aNeighbor_new.append(lNeighbor)
        
        pCell.aNeighbor = aNeighbor_new
        pCell.nNeighbor = len(aNeighbor_new)
        pCell.nNeighbor_land= len(aNeighbor_new)
        pCell.aNeighbor_land = aNeighbor_new
        pCell.nNeighbor_ocean = pCell.nVertex - pCell.nNeighbor_land
        aDggrid_out.append(pCell)

    #calculate neighbor distance
    for pDggrid in aDggrid_out:
        aNeighbor = pDggrid.aNeighbor
        pDggrid.aNeighbor_distance=list()
        for lCellID1 in aNeighbor:
            lIndex = aDggrid_dict[lCellID1]
            pDggrid1 = aDggrid_out[lIndex]  
            dDistance = pDggrid.pVertex_center.calculate_distance( pDggrid1.pVertex_center )
            pDggrid.aNeighbor_distance.append(dDistance)
     
    pDataset = pLayer = pFeature  = None  
    return aDggrid_out

def create_dggrid_mesh(iFlag_global,
                         iFlag_save_mesh,
                         sFilename_mesh,
                         sWorkspace_output, #for dggrid
                         iResolution_index_in = None,
                         sDggrid_type_in= None,
                         iFlag_antarctic_in=None,
                         sFilename_boundary_in = None ):

    #use dggrid table to determine the resolution index   
    sFilename_cell = sWorkspace_output + slash + 'cells'    

    if iResolution_index_in is not None:
        iResolution_index = iResolution_index_in
    else:
        iResolution_index = 10

    if sDggrid_type_in is not None:
        sDggrid_type = sDggrid_type_in
    else:
        sDggrid_type = 'ISEA3H' #default

    dResolution= dggrid_find_resolution_by_index(sDggrid_type, iResolution_index)  
    print('Resolution is: ', dResolution)
    #  
    sResolution = "{:0d}".format( iResolution_index )

    if sFilename_boundary_in is not None:
        if os.path.isfile(sFilename_boundary_in):
            iFlag_crop =1
            sFilename_crop_geojson = sFilename_boundary_in
        else:
            iFlag_crop = 0
    else:
        iFlag_crop = 0   

    iFlag_mode = 1 

    if iFlag_mode == 1: #call the binary directly    
        #write configuration
        sFilename_config= sWorkspace_output + slash +   'dggrid.ini'
        ofs = open(sFilename_config, 'w')
        sLine = 'dggrid_operation GENERATE_GRID' + '\n'
        ofs.write(sLine)
        sLine = 'dggs_type ' + sDggrid_type.upper() + '\n'
        ofs.write(sLine)
        sLine = 'dggs_res_spec ' + sResolution + '\n'
        ofs.write(sLine)

        if iFlag_crop ==1:
            sLine = 'clip_region_files ' + sFilename_crop_geojson + '\n'
            ofs.write(sLine)
            sLine = 'clip_subset_type GDAL'  + '\n'
            ofs.write(sLine)

        else:
            pass

        sLine = 'update_frequency 10000000'+ '\n'
        ofs.write(sLine)
        sLine = 'cell_output_type GDAL_COLLECTION' + '\n'
        ofs.write(sLine)
        sLine = 'cell_output_file_name ' + sFilename_cell + '\n'
        ofs.write(sLine)
        sLine = 'densification 0' + '\n'
        ofs.write(sLine)
        sLine = 'max_cells_per_output_file 0'  + '\n'
        ofs.write(sLine)
        sLine = 'neighbor_output_type GDAL_COLLECTION'  + '\n'
        ofs.write(sLine)       
        
        ofs.close()
        #writen normal run script
        generate_bash_script(sWorkspace_output)
        os.chdir(sWorkspace_output)
        sCommand = "./run_dggrid.sh"
        print(sCommand)
        p = subprocess.Popen(sCommand, shell= True)
        p.wait()

        #convert the pyflowline mesh format
             
        aDggrid = convert_dggrid_mesh_to_pyflowline_mesh(sFilename_cell, sFilename_mesh)
        
    return aDggrid

if __name__ == '__main__':


    sRegion = 'conus'

    sWorkspace_job = '/qfs/people/liao313/jobs/' + 'dggrid' + slash + sRegion + slash + 'simulation'
    sWorkspace_output = '/pic/scratch/liao313/04model/dggrid/' + slash + sRegion + slash + 'simulation'


    sFilename_boundary = '/qfs/people/liao313/data/dggrid/conus/vector/conus_simple.shp'
    iCase_index = 1
    sDggrid_type='isea3h'
    for i  in np.arange(1, 15):
        sWalltime =  "{:02d}".format( iCase_index )
        iResolution_index = i
        create_dggrid_mesh( iCase_index, \
                            iResolution_index ,  \
                            sDggrid_type,

                       sFilename_boundary_in=  sFilename_boundary)






    pass
