import os
import math
import importlib.util
import numpy as np
from osgeo import ogr, osr, gdal
from pyflowline.formats.convert_attributes import convert_gcs_attributes_to_cell
from pyflowline.mesh.jigsaw.run_jigsaw import run_jigsaw #_mpas_workflow import run_jigsaw_mpas_workflow
from pyflowline.algorithms.potentiometric.calculate_potentiometric import calculate_potentiometric
gdal.UseExceptions()
iFlag_cython = importlib.util.find_spec("cython")
from pyearth.gis.geometry.convert_longitude_range import convert_360_to_180_np
from pyearth.gis.geometry.convert_idl_polygon_to_valid_polygon import convert_idl_polygon_to_valid_polygon

def create_mpas_mesh(sFilename_output_in,
        iFlag_global_in = None,
        iFlag_use_mesh_dem_in = None,
        iFlag_save_mesh_in = None,
        iFlag_run_jigsaw_in=None,
        iFlag_antarctic_in=None,
        iFlag_arctic_in=None,
        iFlag_fill_hole_in=None,
        pBoundary_in = None,
        sWorkspace_jigsaw_in = None,
        sFilename_mpas_mesh_netcdf_in= None,
        sFilename_jigsaw_mesh_netcdf_in= None,
        sFilename_land_ocean_mask_in= None,
        aConfig_jigsaw_in = None,
        aFilename_river_network_in = None,
        aFilename_watershed_boundary_in = None,
        aFilename_lake_boundary_in = None,
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
        sFilename_mpas_mesh_netcdf_in (_type_): _description_
        sFilename_output_in (_type_): _description_

    Returns:
        _type_: _description_
    """
    import netCDF4 as nc

    if iFlag_global_in is None:
        iFlag_global=0
    else:
        iFlag_global=iFlag_global_in

    if iFlag_use_mesh_dem_in is None:
        iFlag_use_mesh_dem=0
    else:
        iFlag_use_mesh_dem=iFlag_use_mesh_dem_in

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

    if sFilename_mpas_mesh_netcdf_in is None:
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
    pSpatial_reference_gcs.ImportFromEPSG(4326)  # WGS84 lat/lon
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
        if iFlag_use_mesh_dem == 1:
            pLayer.CreateField(ogr.FieldDefn('elevation_mean', ogr.OFTReal)) #float type for high resolution
            pLayer.CreateField(ogr.FieldDefn('elevation_profile0', ogr.OFTReal)) #float type for high resolution
        else:
            pass

        pLayerDefn = pLayer.GetLayerDefn()
        pFeature = ogr.Feature(pLayerDefn)

    if iFlag_run_jigsaw_in ==1:
        #sFilename_mpas_mesh_netcdf_in = run_jigsaw_mpas_workflow(sWorkspace_jigsaw_in,
        #                         aConfig_in = aConfig_jigsaw_in,
        #    aFilename_river_network_in = aFilename_river_network_in,
        #    aFilename_watershed_boundary_in = aFilename_watershed_boundary_in,
        #    aFilename_lake_boundary_in = aFilename_lake_boundary_in,
        #    aFilename_coastline_in = aFilename_coastline_in)
        projector=[0.0, 0.0]
        geom, gprj, mesh, mprj = run_jigsaw(sWorkspace_jigsaw_in, projector,
                                      aConfig_in=aConfig_jigsaw_in,
                                      aFilename_river_network_in=aFilename_river_network_in,
                                      aFilename_watershed_boundary_in= aFilename_watershed_boundary_in,
                                      aFilename_lake_boundary_in = aFilename_lake_boundary_in,
                                      aFilename_coastline_in = aFilename_coastline_in)

#-------------------------------------- write output for ESM

        iFlag_mpas_tool = 1
        if iFlag_mpas_tool == 1:
            from pyflowline.mesh.jigsaw.saveesm import saveesm
            sFilename_culled_mesh, sFilename_invert_mesh = saveesm(sWorkspace_jigsaw_in, geom, mesh,
                                                                   sFilename_jigsaw_mesh_netcdf_in=sFilename_jigsaw_mesh_netcdf_in,
                                                                   sFilename_land_ocean_mask_in = sFilename_land_ocean_mask_in)
            print('The generated MPAS mesh is: ', sFilename_invert_mesh)
        else:
            #we will a new function to convert jigsaw mesh to mpas mesh11
            print('The algorithm is not completed yet')
            pass

        sFilename_mpas_mesh_netcdf_in = sFilename_invert_mesh

        pass
    else:
        #check whether the jigsaw mesh file exists
        if sFilename_jigsaw_mesh_netcdf_in is not None:
            if os.path.exists(sFilename_jigsaw_mesh_netcdf_in):
                #can we use the jigsaw mesh to create the mpas mesh?
                pass
            else:
                print('JIGSAW Mesh file does not exist!')
        pass

    if (os.path.exists(sFilename_mpas_mesh_netcdf_in)):
        pass
    else:
        print('Mesh file does not exist!')
        return

    pDatasets_in = nc.Dataset(sFilename_mpas_mesh_netcdf_in)
    netcdf_format = pDatasets_in.file_format
    #read new netcdf
    for sKey, aValue in pDatasets_in.variables.items():
        #we need to filter out unused grids based on mpas specs
        if sKey == 'latCell':
            latCell0 = aValue
        else:
            pass
        if sKey == 'lonCell':
            lonCell0 = aValue
        else:
            pass

        if sKey == 'edgesOnCell':
            edgesOnCell0 = aValue
        else:
            pass

        if sKey == 'cellsOnCell':
            cellsOnCell0 = aValue
        else:
            pass

        if sKey == 'cellsOnEdge':
            cellsOnEdge0 = aValue
        else:
            pass

        if sKey == 'verticesOnCell':
            verticesOnCell0 = aValue
        else:
            pass

        if sKey == 'verticesOnEdge':
            verticesOnEdge0 = aValue
        else:
            pass

        if sKey == 'indexToCellID':
            indexToCellID0 = aValue
        else:
            pass

        if sKey == 'indexToEdgeID':
            indexToEdgeID0 = aValue
        else:
            pass

        if sKey == 'indexToVertexID':
            indexToVertexID0 = aValue
        else:
            pass

        if sKey == 'lonVertex':
            lonVertex0 = aValue
        else:
            pass

        if sKey == 'latVertex':
            latVertex0 = aValue
        else:
            pass

        if sKey == 'areaCell':
            areaCell0 = aValue
        else:
            pass

        if sKey == 'bed_elevation':
            bed_elevation0 = aValue
            iFlag_bed_elevation = 1
        else:
            pass

        if sKey == 'ice_thickness':
            ice_thickness0 = aValue
            iFlag_ice_thickness = 1
        else:
            pass

        if sKey == 'areaCell':
            areaCell0 = aValue
        else:
            pass

        if sKey == 'dcEdge':
            dcEdge0 = aValue
        else:
            pass

        if sKey == 'bed_elevation_profile':
            iFlag_elevation_profile = 1
            bed_elevation_profile0 = aValue
        else:
            pass

    aLatitudeVertex = latVertex0[:] / math.pi * 180
    aLongitudeVertex = lonVertex0[:] / math.pi * 180
    #convert unit
    aLatitudeCell = latCell0[:] / math.pi * 180
    aLongitudeCell = lonCell0[:] / math.pi * 180
    aCellsOnCell = cellsOnCell0[:]
    #aCellOnEdge = cellsOnEdge0[:]
    aEdgesOnCell = edgesOnCell0[:]
    aVertexOnCell = verticesOnCell0[:]
    aVertexOnEdge0 = verticesOnEdge0[:]
    aIndexToCellID = indexToCellID0[:]
    #aIndexToEdgeID = indexToEdgeID0[:]
    #aIndexToVertexID = indexToVertexID0[:]
    if iFlag_bed_elevation == 1:
        aBed_elevation = bed_elevation0[:]
    if iFlag_ice_thickness == 1:
        aIce_thickness = ice_thickness0[:]
    aCellArea = areaCell0[:]
    aDcEdge = dcEdge0[:]
    if iFlag_elevation_profile == 1:
        aBed_elevation_profile = bed_elevation_profile0[:]  #elevation
    ncell = len(aIndexToCellID)
    aMpas = list()
    aMpas_dict = dict()
    lCellIndex=0

    aLongitudeCell_180 = convert_360_to_180_np(aLongitudeCell)
    aLongitudeVertex_180 = convert_360_to_180_np(aLongitudeVertex)
    #add a mpas cell into a list
    def add_cell_into_list(aList, i, lCellID, dArea, dElevation_mean, dElevation_profile0, aCoords_gcs  ):
        iFlag_success = 1
        dLongitude_center =  float(aLongitudeCell_180[i])
        dLatitude_center =  float(aLatitudeCell[i])
        if dLongitude_center > 180:
            print('Warning: longitude > 180')
        #vertex
        aCellOnCellIndex = np.array(aCellsOnCell[i,:])
        aEdgesOnCellIndex = np.array(aEdgesOnCell[i,:])
        aVertexOnCellIndex = np.array(aVertexOnCell[i,:])
        dummy0 = np.where(aVertexOnCellIndex > 0)
        aVertexIndex = aVertexOnCellIndex[dummy0]
        if len(aVertexIndex) != len(set(aVertexIndex)):
            print("Duplicates found in aVertexIndex")
            iFlag_success = 0
            return iFlag_success ,aList

        dummy1 = np.where(aEdgesOnCellIndex > 0)
        aEdgeIndex= aEdgesOnCellIndex[dummy1]
        if len(aEdgeIndex) != len(set(aEdgeIndex)):
            print("Duplicates found in aEdgeIndex")
            iFlag_success = 0
            return iFlag_success, aList

        dummy2 = np.where(aCellOnCellIndex > 0)
        aNeighborIndex= (aCellOnCellIndex[dummy2]).astype(int)
        aVertexIndexOnEdge = np.array(aVertexOnEdge0[aEdgeIndex-1,:]).astype((int))

        #check dimensions are consistent
        if len(aVertexIndex) == len(aEdgeIndex) and len(aEdgeIndex) == len(aVertexIndexOnEdge):
            pmpas = convert_gcs_attributes_to_cell(4, dLongitude_center, dLatitude_center, aCoords_gcs, aVertexIndex, aEdgeIndex, aVertexIndexOnEdge)
            if pmpas is None:
                print('Warning: pmpas is None')
                iFlag_success = 0
                return iFlag_success, aList
            pmpas.dArea = dArea
            pmpas.calculate_edge_length()
            pmpas.dLength_flowline = pmpas.dLength_edge #Default
            pmpas.lCellID = lCellID
            pmpas.dElevation_mean  = dElevation_mean
            pmpas.dElevation_profile0 = dElevation_profile0
            #now setup the neighbor information
            pmpas.aNeighbor=aNeighborIndex
            pmpas.nNeighbor=len(aNeighborIndex)
            if pmpas.nNeighbor != pmpas.nVertex:  #this cell is next to the ocean boundary
                pmpas.nNeighbor_land = pmpas.nNeighbor
                pmpas.nNeighbor_ocean = pmpas.nVertex - pmpas.nNeighbor
                pmpas.aNeighbor_land=aNeighborIndex
                pmpas.nNeighbor_land=len(aNeighborIndex)
            else: #this cell is not at the the land-ocean mask coastal line
                pmpas.nNeighbor_land = pmpas.nNeighbor
                pmpas.nNeighbor_ocean = 0
                pmpas.aNeighbor_land = aNeighborIndex
                pmpas.nNeighbor_land=len(aNeighborIndex)

            aDistance=list()
            for j in range(pmpas.nNeighbor):
                #find shared edge
                lEdgeID= aEdgeIndex[j]
                lIndex = lEdgeID-1
                dDistance = aDcEdge[lIndex]
                aDistance.append(dDistance)
                pass

            #this contains all the original mpas neighbor distance
            pmpas.aNeighbor_distance = aDistance
            aList.append(pmpas)
            return iFlag_success, aList
        else:
            print('Warning: len(aVertexIndex) != len(aVertexIndexOnEdge)', 'cellID:', lCellID)
            #we will not add this cell if something is wrong
            iFlag_success = 0
            return iFlag_success, aList

    if iFlag_antarctic == 1: #use potentiometric
        #if it is antarctic, we dont need the boundary
        for i in range(ncell):
            #center
            dLongitude_center = float(aLongitudeCell_180[i])
            dLatitude_center = float(aLatitudeCell[i])
            aVertexOnCellIndex = np.array(aVertexOnCell[i,:])
            dummy0 = np.where(aVertexOnCellIndex > 0)
            aVertexIndex = aVertexOnCellIndex[dummy0]
            aLonVertex = aLongitudeVertex_180[aVertexIndex-1]
            aLatVertex = aLatitudeVertex[aVertexIndex-1]
            nVertex = len(aLonVertex)
            if nVertex < 3:
                print('Warning: nVertex < 3')
                continue
            #first check if it is within the boundary
            iFlag = False
            ring = ogr.Geometry(ogr.wkbLinearRing)
            aCoords_gcs = np.full((nVertex,2), -9999.0, dtype=float)
            for j in range(nVertex):
                x1 = aLonVertex[j]
                y1 = aLatVertex[j]
                ring.AddPoint(x1, y1)
                aCoords_gcs[j,0] = x1
                aCoords_gcs[j,1] = y1
                pass

            x1 = aLonVertex[0]
            y1 = aLatVertex[0]
            ring.AddPoint(x1, y1) #double check
            pPolygon = ogr.Geometry(ogr.wkbPolygon)
            pPolygon.AddGeometry(ring)

            if dLatitude_center < -60:
                iFlag = True
            else:
                iFlag = False
                pass

            if ( iFlag == True ):
                #check polygon
                if pPolygon.IsValid() == False:
                    print('Warning: invalid polygon')
                    continue
                lCellID = int(aIndexToCellID[i])
                if iFlag_bed_elevation == 1:
                    dElevation_mean = float(aBed_elevation[i])
                else:
                    dElevation_mean = 0.0 #-9999
                if iFlag_elevation_profile == 1:
                    dElevation_profile0 = float(aBed_elevation_profile[i,0])
                else:
                    dElevation_profile0 = 0.0 #-9999
                if iFlag_ice_thickness == 1:
                    dThickness_ice = float( aIce_thickness[i] )
                else:
                    dThickness_ice = 0.0 #-9999
                dArea = float(aCellArea[i])

                #then check if it is ice free

                if dThickness_ice > 0 :
                    #use potentiometric
                    dElevation_mean = calculate_potentiometric(dElevation_mean , dThickness_ice)
                    dElevation_profile0 = calculate_potentiometric(dElevation_profile0 , dThickness_ice)

                #call fuction to add the cell
                iFlag_success, aMpas = add_cell_into_list(aMpas, i, lCellID, dArea, dElevation_mean, dElevation_profile0, aCoords_gcs )
                if iFlag_success == 1:
                    aMpas_dict[lCellID] = lCellIndex
                    lCellIndex = lCellIndex + 1
                #save mesh cell
                if iFlag_save_mesh ==1:
                    pFeature.SetGeometry(pPolygon)
                    pFeature.SetField("cellid", int(lCellID) )
                    pFeature.SetField("longitude", dLongitude_center )
                    pFeature.SetField("latitude", dLatitude_center )
                    pFeature.SetField("area", dArea )
                    if iFlag_use_mesh_dem == 1:
                        pFeature.SetField("elevation_mean", dElevation_mean )
                        pFeature.SetField("elevation_profile0", dElevation_profile0 )

                    pLayer.CreateFeature(pFeature)
    else:
        if iFlag_arctic == 1: #for arctic only
            for i in range(ncell):
                #center
                dLongitude_center = float(aLongitudeCell_180[i])
                dLatitude_center = float(aLatitudeCell[i])
                aVertexOnCellIndex = np.array(aVertexOnCell[i,:])
                dummy0 = np.where(aVertexOnCellIndex > 0)
                aVertexIndex = aVertexOnCellIndex[dummy0]
                aLonVertex = aLongitudeVertex_180[aVertexIndex-1]
                aLatVertex = aLatitudeVertex[aVertexIndex-1]
                nVertex = len(aLonVertex)
                if nVertex < 3:
                    print('Warning: nVertex < 3')
                    continue
                #first check if it is within the boundary
                iFlag = False
                ring = ogr.Geometry(ogr.wkbLinearRing)
                aCoords_gcs = np.full((nVertex,2), -9999.0, dtype=float)
                for j in range(nVertex):
                    x1 = aLonVertex[j]
                    y1 = aLatVertex[j]
                    ring.AddPoint(x1, y1)
                    aCoords_gcs[j,0] = x1
                    aCoords_gcs[j,1] = y1
                    pass

                x1 = aLonVertex[0]
                y1 = aLatVertex[0]
                ring.AddPoint(x1, y1) #double check
                pPolygon = ogr.Geometry(ogr.wkbPolygon)
                pPolygon.AddGeometry(ring)

                if dLatitude_center > 55.0:
                    iFlag = True
                else:
                    iFlag = False
                    pass

                if ( iFlag == True ):
                    if pPolygon.IsValid() == False:
                        print('Warning: invalid polygon')
                        continue
                    lCellID = int(aIndexToCellID[i])
                    if iFlag_bed_elevation == 1:
                        dElevation_mean = float(aBed_elevation[i])
                    else:
                        dElevation_mean = 0.0 #-9999
                    if iFlag_elevation_profile == 1:
                        dElevation_profile0 = float(aBed_elevation_profile[i,0])
                    else:
                        dElevation_profile0 = 0.0 #-9999
                    if iFlag_ice_thickness == 1:
                        dThickness_ice = float( aIce_thickness[i] )
                    else:
                        dThickness_ice = 0.0 #-9999
                    dArea = float(aCellArea[i])

                    #then check if it is ice free
                    if dThickness_ice > 0 :
                        #use potentiometric
                        dElevation_mean = calculate_potentiometric(dElevation_mean , dThickness_ice)
                        dElevation_profile0 = calculate_potentiometric(dElevation_profile0 , dThickness_ice)


                    #call fuction to add the cell
                    iFlag_success, aMpas = add_cell_into_list(aMpas, i, lCellID, dArea, dElevation_mean, dElevation_profile0, aCoords_gcs )
                    if iFlag_success == 1:
                        aMpas_dict[lCellID] = lCellIndex
                        lCellIndex = lCellIndex + 1
                    #save mesh cell
                    if iFlag_save_mesh ==1:
                        pFeature.SetGeometry(pPolygon)
                        pFeature.SetField("cellid", int(lCellID) )
                        pFeature.SetField("longitude", dLongitude_center )
                        pFeature.SetField("latitude", dLatitude_center )
                        pFeature.SetField("area", dArea )
                        if iFlag_use_mesh_dem == 1:
                            pFeature.SetField("elevation_mean", dElevation_mean )
                            pFeature.SetField("elevation_profile0", dElevation_profile0 )

                        pLayer.CreateFeature(pFeature)
            pass
        else:
            for i in range(ncell):
                dLongitude_center = float(aLongitudeCell_180[i])
                dLatitude_center = float(aLatitudeCell[i])
                #vertex
                aVertexOnCellIndex = np.array(aVertexOnCell[i,:])
                dummy0 = np.where(aVertexOnCellIndex > 0)
                aVertexIndex = aVertexOnCellIndex[dummy0]
                aLonVertex = aLongitudeVertex_180[aVertexIndex-1]
                aLatVertex = aLatitudeVertex[aVertexIndex-1]
                nVertex = len(aLonVertex)
                if nVertex < 3:
                    print('Warning: nVertex < 3')
                    continue
                #first check if it is within the boundary
                ring = ogr.Geometry(ogr.wkbLinearRing)
                aCoords_gcs = np.full((nVertex,2), -9999.0, dtype=float)
                for j in range(nVertex):
                    x1 = aLonVertex[j]
                    y1 = aLatVertex[j]
                    ring.AddPoint(x1, y1)
                    aCoords_gcs[j,0] = x1
                    aCoords_gcs[j,1] = y1
                    pass

                x1 = aLonVertex[0]
                y1 = aLatVertex[0]
                ring.AddPoint(x1, y1) #double check
                pPolygon = ogr.Geometry(ogr.wkbPolygon)
                pPolygon.AddGeometry(ring)

                #check within first
                dLon_min = np.min(aCoords_gcs[:,0])
                dLon_max = np.max(aCoords_gcs[:,0])
                iFlag_debug = 0
                #if iFlag_debug == 1:
                #    if dLongitude_center > -90 and dLongitude_center < -80 and dLatitude_center > 40 and dLatitude_center < 50:
                #        print('test mesh cell')
                #    else:
                #        pass
                if iFlag_global == 1:
                    if dLatitude_center >= -60:
                        if dLatitude_center>85: #exclude the north pole as well
                            iFlag = False
                            continue
                        else:
                            if np.abs(dLon_min-dLon_max) > 100: #this polygon cross international date line
                                #print('Warning: longitude > 180')
                                #keep this but use an alternative way to check
                                iFlag = True #do not check the polygon?
                                pPolygon_new = convert_idl_polygon_to_valid_polygon(pPolygon)
                                if pPolygon_new is not None:
                                    if pPolygon_new.IsValid() == False:
                                        print('Warning: invalid polygon')
                                        continue
                                    else:
                                        pass
                                else:
                                    continue

                            else:
                                iFlag = True
                                if pPolygon.IsValid() == False:
                                    print('Warning: invalid polygon')
                                    continue
                    else:
                        iFlag = False  #remove antiarctic from global mesh for now
                        continue
                else:
                    iFlag = False
                    if np.abs(dLon_min-dLon_max) > 100: #this polygon cross international date line
                        #print('Warning: longitude > 180')
                        continue
                    else:
                        if pPolygon.IsValid() == False:
                            print('Warning: invalid polygon')
                            continue
                        else:
                            if pPolygon.Within(pBoundary):
                                iFlag = True
                            else:
                                #then check intersection
                                if pPolygon.Intersects(pBoundary):
                                    iFlag = True


                if ( iFlag == True ):
                    lCellID = int(aIndexToCellID[i])
                    if iFlag_bed_elevation == 1:
                        dElevation_mean = float(aBed_elevation[i])
                    else:
                        dElevation_mean = 0.0 #-9999
                    if iFlag_elevation_profile == 1:
                        dElevation_profile0 = float(aBed_elevation_profile[i,0])
                    else:
                        dElevation_profile0 = 0.0
                    if iFlag_ice_thickness == 1:
                        dThickness_ice = float( aIce_thickness[i] )
                    else:
                        dThickness_ice = 0.0 #-9999

                    #then check if it is ice free
                    if dThickness_ice > 0 :
                        #use potentiometric
                        dElevation_mean = calculate_potentiometric(dElevation_mean , dThickness_ice)
                        if iFlag_elevation_profile == 1:
                            dElevation_profile0 = calculate_potentiometric(dElevation_profile0 , dThickness_ice)

                    dArea = float(aCellArea[i])
                    #then check if it is ice free
                    #if iFlag_remove_ice == 1:
                    #    if dThickness_ice > 0 :
                    #        continue
                    #    else:
                    #        pass
                    #else:
                    #    pass
                    #call fuction to add the cell

                    iFlag_success, aMpas = add_cell_into_list(aMpas, i, lCellID, dArea, dElevation_mean, dElevation_profile0, aCoords_gcs )
                    if iFlag_success == 1:
                        aMpas_dict[lCellID] = lCellIndex
                        lCellIndex = lCellIndex + 1
                    #save mesh cell
                    if iFlag_save_mesh == 1:
                        dLongitude_center = float(aLongitudeCell_180[i])
                        dLatitude_center = float(aLatitudeCell[i])
                        pFeature.SetGeometry(pPolygon)
                        pFeature.SetField("cellid", int(lCellID) )
                        pFeature.SetField("longitude", dLongitude_center )
                        pFeature.SetField("latitude", dLatitude_center )
                        pFeature.SetField("area", dArea )
                        if iFlag_use_mesh_dem == 1:
                            pFeature.SetField("elevation_mean", dElevation_mean )
                            pFeature.SetField("elevation_profile0", dElevation_profile0 )

                        pLayer.CreateFeature(pFeature)

    #close the dataset
    if iFlag_save_mesh == 1:
        pFeature = None
        pLayer = None
        pDataset = None

    #for maps we need to clean some cell because they were not actually in the domain
    #besides, we need to add some smal holes back
    #to do this, we need two steps.
    if iFlag_global == 1:
        aMpas_out = aMpas
    else:
        aMpas_out = list()
        ncell = len(aMpas)
        #still need to get rid cell that are not in the domain
        for pCell in aMpas:
            aNeighbor = pCell.aNeighbor
            aNeighbor_land_update = list()
            for lNeighbor in aNeighbor:
                if lNeighbor in aMpas_dict:
                    aNeighbor_land_update.append(lNeighbor)
            #for latlon, there is no ocean concept
            pCell.aNeighbor = aNeighbor_land_update
            pCell.nNeighbor= len(aNeighbor_land_update)
            pCell.aNeighbor_land = aNeighbor_land_update
            pCell.nNeighbor_land= len(aNeighbor_land_update)
            pCell.nNeighbor_ocean = pCell.nVertex - pCell.nNeighbor_land
            aMpas_out.append(pCell)

    return aMpas_out


