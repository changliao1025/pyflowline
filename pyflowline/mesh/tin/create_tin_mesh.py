import os
import math
import importlib.util
import numpy as np
from osgeo import ogr, osr, gdal
from pyearth.gis.location.xyz_to_lonlat import xyz_to_lonlat
from pyflowline.formats.convert_coordinates import convert_gcs_coordinates_to_cell
from pyflowline.formats.convert_attributes import convert_gcs_attributes_to_cell
from pyflowline.mesh.jigsaw.run_jigsaw import run_jigsaw
from pyflowline.mesh.jigsaw.savetin import savetin
from pyflowline.algorithms.potentiometric.calculate_potentiometric import calculate_potentiometric
gdal.UseExceptions()
iFlag_cython = importlib.util.find_spec("cython")
if iFlag_cython is not None:
    from pyflowline.algorithms.cython.kernel import convert_360_to_180
else:
    from pyearth.gis.geometry.convert_longitude_range import convert_360_to_180

def update_vertices_on_cell(nCell, nVertex, cellsOnVertex):
    # Initialize verticesOnCell with empty lists for each cell
    verticesOnCell = [[] for _ in range(nCell)]
    verticesOnCellSets = [set() for _ in range(nCell)]  # Use sets for faster lookup

    # Iterate over each vertex
    for iVertex in range(nVertex):
        # Iterate over each cell connected to the current vertex
        for iCell in cellsOnVertex[iVertex]:
            if iCell != -1:
                cellIndex = iCell - 1
                # Check if the vertex is already in the set for the cell
                if iVertex not in verticesOnCellSets[cellIndex]:
                    verticesOnCell[cellIndex].append(iVertex)
                    verticesOnCellSets[cellIndex].add(iVertex)

    return verticesOnCell

def create_tin_mesh(   sFilename_output_in,
                    iFlag_global_in = None,
    iFlag_save_mesh_in= None,
    sFilename_jigsaw_mesh_netcdf_in=None,
    iFlag_run_jigsaw_in=None,
    iFlag_antarctic_in=None,
    iFlag_arctic_in=None,
    pBoundary_in = None,
    sWorkspace_jigsaw_in = None,
    aConfig_jigsaw_in = None,
    aFilename_river_in = None,
    aFilename_watershed_boundary_in = None,
    aFilenamae_lake_boundary_in = None,
    aFilename_coastline_in = None,
    iFlag_read_mesh_in  = None,
    iFlag_generate_mesh_in=None):
    """
    Create a MPAS mesh

    Args:
        iFlag_global_in (int): _description_
        iFlag_use_mesh_dem (int): _description_
        iFlag_save_mesh_in (int): _description_
        pBoundary_in (_type_): _description_
        sFilename_jigsaw_mesh_netcdf_in (_type_): _description_
        sFilename_output_in (_type_): _description_

    Returns:
        _type_: _description_
    """
    import netCDF4 as nc

    if iFlag_global_in is None:
        iFlag_global=0
    else:
        iFlag_global=iFlag_global_in


    if iFlag_save_mesh_in is None:
        iFlag_save_mesh=1
    else:
        iFlag_save_mesh=iFlag_save_mesh_in

    if iFlag_antarctic_in is None:
        iFlag_antarctic=0
    else:
        iFlag_antarctic=iFlag_antarctic_in

    if iFlag_arctic_in is None:
        iFlag_arctic=0
    else:
        iFlag_arctic=iFlag_arctic_in

    if sFilename_jigsaw_mesh_netcdf_in is None:
        print('This mesh file will be generated!')
    else:
        pass


    if pBoundary_in is None:
        pBoundary = None
    else:
        #for the reason that a geometry object will be crash if the associated dataset is closed, we must pass wkt string
        #https://gdal.org/api/python_gotchas.html
        pBoundary = ogr.CreateGeometryFromWkt(pBoundary_in)
        #pBoundary = pBoundary_in #only works for shapely geometry object

    if os.path.exists(sFilename_output_in):
        os.remove(sFilename_output_in)

    iFlag_elevation_profile=0
    iFlag_bed_elevation=0
    iFlag_ice_thickness=0
    pDriver_geojson = ogr.GetDriverByName('GeoJSON')
    pSpatial_reference_gcs = osr.SpatialReference()
    pSpatial_reference_gcs.ImportFromEPSG(4326)    # WGS84 lat/lon
    #geojson
    if iFlag_save_mesh ==1:
        pDataset = pDriver_geojson.CreateDataSource(sFilename_output_in)
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

    if iFlag_run_jigsaw_in ==1:
        projector=[0.0, 0.0]
        geom, gprj, mesh, mprj =  run_jigsaw(sWorkspace_jigsaw_in, projector,
                                      aConfig_in=aConfig_jigsaw_in,
                                      aFilename_river_in=aFilename_river_in,
                                      aFilename_watershed_boundary_in= aFilename_watershed_boundary_in,
                                      aFilenamae_lake_boundary_in = aFilenamae_lake_boundary_in,
                                      aFilename_coastline_in = aFilename_coastline_in)

#-------------------------------------- write output for ESM

        sFilename_jigsaw_mesh_netcdf_in = savetin(sWorkspace_jigsaw_in, geom, mesh)

        pass
    else:
        pass

    if (os.path.exists(sFilename_jigsaw_mesh_netcdf_in)):
        pass
    else:
        print('Mesh file does not exist!')
        return

    pDatasets_in = nc.Dataset(sFilename_jigsaw_mesh_netcdf_in)

    netcdf_format = pDatasets_in.file_format
    #read new netcdf
    for sKey, aValue in pDatasets_in.variables.items():
        #we need to filter out unused grids based on mpas specs
        if sKey == 'cellsOnVertex':
            aCellOnVertex = aValue
        else:
            pass
        if sKey == 'xVertex':
            aX = aValue
        else:
            pass

        if sKey == 'yVertex':
            aY = aValue
        else:
            pass

        if sKey == 'zVertex':
            aZ = aValue
        else:
            pass

        if sKey == 'xCell':
            aX0 = aValue
        else:
            pass
        if sKey == 'yCell':
            aY0 = aValue
        else:
            pass
        if sKey == 'zCell':
            aZ0 = aValue
        else:
            pass

    #build connectivity
    nCell = len(aX0)
    nVertex = len(aX)
    #verticesOnCell = update_vertices_on_cell(nCell, nVertex, aCellOnVertex)

    aTin = list()
    aTin_dict = dict()
    lCellIndex=0

    if iFlag_antarctic == 1: #use potentiometric
        iFlag_remove_ice = 0
        #if it is antarctic, we dont need the boundary
        pass
    else:
        if iFlag_arctic ==1: #for arctic only
            iFlag_remove_ice = 0
            pass
        else:
            lCellID = 1
            for i in range(nVertex):
                dummy = aCellOnVertex[i]
                nVertex1 = len(dummy)
                ring = ogr.Geometry(ogr.wkbLinearRing)
                aCoords_gcs = np.full((nVertex1+1,2), -9999.0, dtype=float)
                for j in range(nVertex1):
                    dummy1 = dummy[j]
                    x, y, z = aX0[dummy1 - 1], aY0[dummy1 - 1], aZ0[dummy1 - 1]
                    x1, y1 = xyz_to_lonlat(x, y, z)
                    ring.AddPoint(x1, y1)
                    aCoords_gcs[j] = [x1, y1]

                # Close the polygon by adding the first point again
                x, y, z = aX0[dummy[0] - 1], aY0[dummy[0] - 1], aZ0[dummy[0] - 1]
                x1, y1 = xyz_to_lonlat(x, y, z)
                aCoords_gcs[nVertex1] = [x1, y1]
                ring.AddPoint(x1, y1)

                pPolygon = ogr.Geometry(ogr.wkbPolygon)
                pPolygon.AddGeometry(ring)

                if iFlag_global == 1:
                    iFlag = True
                else:
                    iFlag = False
                    if pPolygon.Within(pBoundary):
                        iFlag = True
                    else:
                        dLon_min, dLon_max = np.min(aCoords_gcs[:, 0]), np.max(aCoords_gcs[:, 0])
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
                    dLongitude_center =float(np.mean(aCoords_gcs[0:nVertex,0]) ) #these are not the actual center, see circumcenter function
                    dLatitude_center = float(np.mean(aCoords_gcs[0:nVertex,1]))
                    pTin = convert_gcs_coordinates_to_cell(7, dLongitude_center, dLatitude_center, aCoords_gcs)
                    pTin.lCellID = lCellID
                    dArea = pTin.calculate_cell_area()
                    pTin.dArea = dArea
                    pTin.calculate_edge_length()
                    aTin.append(pTin)
                    pFeature.SetGeometry(pPolygon)
                    pFeature.SetField("cellid", lCellID)
                    pFeature.SetField("longitude", dLongitude_center )
                    pFeature.SetField("latitude", dLatitude_center )
                    pFeature.SetField("area", dArea )
                    pLayer.CreateFeature(pFeature)
                    aTin_dict[lCellID] = lCellIndex
                    print('cell id: ', lCellID)
                    lCellIndex = lCellIndex + 1
                    lCellID = lCellID + 1

    aTin_out = aTin
    return aTin_out


