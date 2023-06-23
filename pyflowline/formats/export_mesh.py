import os
import json



def convert_mesh_to_kml(sFilename_geojson_in, sFilename_kml_out):
    """
    Convert the mesh to a kml file
    """
    import simplekml
    with open(sFilename_geojson_in) as f:
        data = json.load(f)
        kml = simplekml.Kml()
        for feature in data['features']:
            geom = feature['geometry']
            geom_type = geom['type']
            lCellID = feature['properties']['lCellID']
            sCellID = str(lCellID)

            if geom_type == 'Polygon':
                kml.newpolygon(name=sCellID,
                               description='cell',
                               outerboundaryis=geom['coordinates'][0])
            elif geom_type == 'LineString':
                kml.newlinestring(name=sCellID,
                                  description='flowline',
                                  coords=geom['coordinates'])
            elif geom_type == 'Point':
                kml.newpoint(name=sCellID,
                             description='point',
                             coords=[geom['coordinates']])
            else:
                print("ERROR: unknown type:", geom_type)
        kml.save(sFilename_kml_out)

    return


def export_mesh_info_to_json(aCell_in, aFlowline_in, aCellID_outlet_iin, sFilename_json_in):
    """
    Export the mesh information into a json file

    Args:
        aCell_in (_type_): _description_
        aFlowline_in (_type_): _description_
        aCellID_outlet_iin (_type_): _description_
        sFilename_json_in (_type_): _description_
    """
    if os.path.exists(sFilename_json_in): 
        os.remove(sFilename_json_in)
        pass
    ncell=len(aCell_in)
    nFlowline = len(aFlowline_in)
    for i in range(nFlowline):
        pFlowline = aFlowline_in[i]
        nEdge = pFlowline.nEdge
        nVertex = pFlowline.nVertex
        aEdge = pFlowline.aEdge
        iStream_segment = pFlowline.iStream_segment
        iStream_order = pFlowline.iStream_order
        for j in range(nEdge):
            pEdge = aEdge[j]
            pVertex_start = pEdge.pVertex_start
            pVertex_end = pEdge.pVertex_end
            for k in range(ncell):
                pVertex_center = aCell_in[k].pVertex_center
                if pVertex_center == pVertex_start:
                    aCell_in[k].iStream_segment_burned = iStream_segment
                    aCell_in[k].iStream_order_burned = iStream_order
                    for l in range(ncell):
                        pVertex_center2 = aCell_in[l].pVertex_center
                        lCellID = aCell_in[l].lCellID
                        if pVertex_center2 == pVertex_end:
                            aCell_in[k].lCellID_downstream_burned = lCellID
                            if lCellID in aCellID_outlet_iin:
                                aCell_in[l].iStream_segment_burned = iStream_segment
                                aCell_in[l].iStream_order_burned = iStream_order
                            break
                        
                        
        pass
    

    with open(sFilename_json_in, 'w', encoding='utf-8') as f:
        sJson = json.dumps([json.loads(ob.tojson()) for ob in aCell_in], indent = 4)        
        f.write(sJson)    
        f.close()
    
    return
