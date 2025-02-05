import os, stat
import numpy as np
from osgeo import ogr, osr
import subprocess
import netCDF4 as nc
from pyearth.system.define_global_variables import *
from pyflowline.classes.vertex import pyvertex
from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_cell

def cube_to_sphere(x, y, z):
    # Normalize the coordinates
    norm = np.sqrt(x**2 + y**2 + z**2)
    x /= norm
    y /= norm
    z /= norm

    # Convert to spherical coordinates
    lon = np.arctan2(y, x)
    lat = np.arcsin(z)

    # Convert radians to degrees
    lon = np.degrees(lon)
    lat = np.degrees(lat)

    return float(lon), float(lat)

def generate_bash_script(sWorkspace_output, sFilename_mesh_in, sResolution_in):

    #example call
    #args = ['GenerateCSMesh', '--alt', '--res', f'{self.resolution}',
    #    '--file', f'{self.resolution_name}.g']
    #subprocess.run(args)

    #sFilename_configuration  =  os.path.join(sWorkspace_output,  sName )

    sCommand = '--alt --res ' + sResolution_in + ' --file ' + sFilename_mesh_in

    os.chdir(sWorkspace_output)
    #detemine the system platform
    # Determine the appropriate executable name for the platform
    system = platform.system()
    if platform.system() == 'Windows':
        sFilename_executable = 'GenerateCSMesh.exe'
        iFlag_unix = 0
    else:
        sFilename_executable = 'GenerateCSMesh'

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
        sFilename_bash = os.path.join(str(Path(  sWorkspace_output)  ) ,  "run_csmesh.sh" )
        ofs = open(sFilename_bash, 'w')
        sLine = '#!/bin/bash\n'
        ofs.write(sLine)
        sLine = 'cd ' +   sWorkspace_output+ '\n'
        ofs.write(sLine)
        sLine = sFilename_executable + ' ' + sCommand + '\n'
        ofs.write(sLine)
        ofs.close()
        os.chmod(sFilename_bash, stat.S_IRWXU )
    else:
        sFilename_bash = os.path.join(str(Path(  sWorkspace_output)  ) ,  "run_dggrid.bat" )
        ofs = open(sFilename_bash, 'w')
        sLine = 'cd ' +   sWorkspace_output+ '\n'
        ofs.write(sLine)
        sLine = sFilename_executable + ' ' + sCommand + '\n'
        ofs.write(sLine)
        ofs.close()
        os.chmod(sFilename_bash, stat.S_IRWXU )

    return

def convert_cubicspheremesh_to_pyflowline_mesh(iFlag_global_in,
                                               sFilename_csmesh_in,
                                               sFilename_mesh_pyflowline,
                                         pBoundary_in=None      ):

    if pBoundary_in is None:
        pBoundary = None
    else:
        #for the reason that a geometry object will be crash if the associated dataset is closed, we must pass wkt string
        #https://gdal.org/api/python_gotchas.html
        pBoundary = ogr.CreateGeometryFromWkt(pBoundary_in)
        #pBoundary = pBoundary_in #only works for shapely geometry object

    pDateset = nc.Dataset(sFilename_csmesh_in, 'r')

    # Print the dimensions
    print("Dimensions:")
    for dim in pDateset.dimensions.values():
        print(dim)

    # Print the variables
    print("\nVariables:")
    for var in pDateset.variables.values():
        print(var)

    # Access specific data (example: coordinates)
    if 'coord' in pDateset.variables:
        coords = pDateset.variables['coord'][:]
        print("\nCoordinates:")

    if 'connect1' in pDateset.variables:
        aConnection = pDateset.variables['connect1'][:]
        print("\nconnect1:")

    if 'global_id1' in pDateset.variables:
        aGlobalID = pDateset.variables['global_id1'][:]
        print("\nglobal_id1:")


    nCell = len(aConnection)

    aCubicsphere = list()
    pDriver_geojson = ogr.GetDriverByName('GeoJSON')
    pDataset = pDriver_geojson.CreateDataSource(sFilename_mesh_pyflowline)
    pSpatial_reference_gcs = osr.SpatialReference()
    pSpatial_reference_gcs.ImportFromEPSG(4326)    # WGS84 lat/lon
    pSpatial_reference_gcs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
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

    lCellIndex = 0
    for i in range(nCell):

        lCellID = aGlobalID[i]

        pConnection = aConnection[i]
        nVertex = len(pConnection)
        ring = ogr.Geometry(ogr.wkbLinearRing)
        aCoords_gcs  = np.full((nVertex + 1,2), -9999.0, dtype=float)
        for j in range(nVertex):
            x = coords[0][pConnection[j]-1]
            y = coords[1][pConnection[j]-1]
            z = coords[2][pConnection[j]-1]
            dLongitude, dLatitude = cube_to_sphere(x, y, z)
            aCoords_gcs[j][0] = dLongitude
            aCoords_gcs[j][1] = dLatitude
            ring.AddPoint(dLongitude, dLatitude)
            pass

        aCoords_gcs[nVertex][0] = aCoords_gcs[0][0]
        aCoords_gcs[nVertex][1] = aCoords_gcs[0][1]
        dLongitude_center = np.mean(aCoords_gcs[0:4,0])
        dLatitude_center = np.mean(aCoords_gcs[0:4,1])

        #add the first vertex to the end to close the polygon
        ring.AddPoint(aCoords_gcs[0][0], aCoords_gcs[0][1])
        pPolygon = ogr.Geometry(ogr.wkbPolygon)
        pPolygon.AddGeometry(ring)
        #check within first
        iFlag = False

        if iFlag_global_in == 1:
            iFlag = True
        else:
            if pPolygon.Within(pBoundary):
                iFlag = True
            else:
                dLon_min = np.min(aCoords_gcs[:,0])
                dLon_max = np.max(aCoords_gcs[:,0])
                if np.abs(dLon_min-dLon_max) > 100: #this polygon cross international date line
                    #print('Warning: longitude > 180')
                    pass
                else:
                    #then check intersection
                    if pPolygon.Intersects(pBoundary):
                        iFlag = True
                    else:
                        pass

        if ( iFlag == True ):
            pcubicsphere = convert_gcs_coordinates_to_cell(8, dLongitude_center, dLatitude_center, aCoords_gcs)
            dArea = pcubicsphere.calculate_cell_area()
            pcubicsphere.calculate_edge_length()
            pcubicsphere.dLength_flowline = pcubicsphere.dLength #Default
            pcubicsphere.lCellID = lCellID
            aCubicsphere.append(pcubicsphere)
            lCellIndex = lCellIndex + 1
            nVertex = pcubicsphere.nVertex
            ring = ogr.Geometry(ogr.wkbLinearRing)
            for j in range(nVertex):
                x1 = pcubicsphere.aVertex[j].dLongitude_degree
                y1 = pcubicsphere.aVertex[j].dLatitude_degree
                ring.AddPoint(x1, y1)
                pass
            x1 = pcubicsphere.aVertex[0].dLongitude_degree
            y1 = pcubicsphere.aVertex[0].dLatitude_degree
            ring.AddPoint(x1, y1) #double check
            pPolygon = ogr.Geometry(ogr.wkbPolygon)
            pPolygon.AddGeometry(ring)
            pFeature.SetGeometry(pPolygon)
            pFeature.SetField("cellid", int(lCellID) )
            pFeature.SetField("longitude", dLongitude_center )
            pFeature.SetField("latitude", dLatitude_center )
            pFeature.SetField("area", dArea )
            pLayer.CreateFeature(pFeature)


    pLayer = pFeature = pDataset = None

    return aCubicsphere

def create_cubicsphere_mesh(iFlag_global_in,
                            iFlag_save_mesh_in,
                       dResolution_meter_in,
                       sFilename_output_in,
                       sWorkspace_output,
                       pBoundary_in=None):
    if pBoundary_in is None:
        print('Are you sure you want to create a mesh without boundary?')
        print('In this case, the program will generate a boundary for you using the bounding box of the mesh')
        #create a boundary for the mesh

    else:
        pBoundary = ogr.CreateGeometryFromWkt(pBoundary_in)

    if os.path.exists(sFilename_output_in):
        os.remove(sFilename_output_in)

    #convert the resolution to km
    dResolution_km = dResolution_meter_in / 1000.0
    #radius of the earth
    dRadius = 6371.0
    #calculate the circumference
    dCircumference = 2.0 * np.pi * dRadius
    dLength = dCircumference / 4.0
    #calculate the number of cells
    nResulution = int(dLength / dResolution_km)
    #convert the resolution to string
    sResolution = str(nResulution)

    print('The resolution is ' + sResolution )

    sFilename_csmesh = sWorkspace_output + slash + 'cubic_sphere_mesh.g'
    generate_bash_script(sWorkspace_output, sFilename_csmesh, sResolution)
    os.chdir(sWorkspace_output)
    sCommand = "./run_csmesh.sh"
    p = subprocess.Popen(sCommand, shell= True)
    p.wait()

    aCubicsphere = convert_cubicspheremesh_to_pyflowline_mesh(iFlag_global_in,
                                                              sFilename_csmesh,
                                                               sFilename_output_in,
                                                                 pBoundary_in = pBoundary_in)

    return aCubicsphere

