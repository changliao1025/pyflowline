from osgeo import ogr
import simplekml
def convert_geojson_to_kml(sFilename_geojson, sFilename_kml):
    # Read the GeoJSON file
    
    #geojson_driver = ogr.GetDriverByName('GeoJSON')
    #geojson_data_source = geojson_driver.Open(sFilename_geojson)
    #layer = geojson_data_source.GetLayer()

    # Create a KML object
    #kml = simplekml.Kml()

    # Iterate over features in the GeoJSON layer
    #feature = layer.GetNextFeature()
    #while feature:
        # Get the geometry of the feature
    ##    geometry = feature.GetGeometryRef()
    #    sGeometry_type = geometry.GetGeometryName()   
        # Convert the geometry to WKT
    #    wkt = geometry.ExportToWkt()
    #    if(sGeometry_type == 'POLYGON'):
            # Create a KML placemark and add it to the KML object
            #polygon = kml.newpolygon()
            #polygon.geometry = simplekml.Geometry(wkt=wkt)
    #        pass
    #    else:
    #        if (sGeometry_type == 'LINESTRING'):
                # Create a KML placemark and add it to the KML object
                #linestring = kml.newlinestring()
                #linestring.geometry = simplekml.Geometry(wkt=wkt)
    #            pass

        # Move to the next feature
        #feature = layer.GetNextFeature()

    # Save the KML to a file
    #sFilename_kml = 'output.kml'
    #kml.save(sFilename_kml)
    return