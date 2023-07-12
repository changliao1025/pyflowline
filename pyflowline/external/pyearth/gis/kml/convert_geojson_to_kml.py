from osgeo import ogr
import simplekml
def convert_geojson_to_kml(sFilename_geojson, sFilename_kml):
    # Read the GeoJSON file
    
    geojson_driver = ogr.GetDriverByName('GeoJSON')
    geojson_data_source = geojson_driver.Open(sFilename_geojson)
    layer = geojson_data_source.GetLayer()

    # Create a KML object
    kml = simplekml.Kml()

    # Iterate over features in the GeoJSON layer
    feature = layer.GetNextFeature()
    while feature:
        # Get the geometry of the feature
        geometry = feature.GetGeometryRef()

        # Convert the geometry to WKT
        wkt = geometry.ExportToWkt()

        # Create a KML placemark and add it to the KML object
        polygon = kml.newpolygon()
        polygon.geometry = simplekml.Geometry(wkt=wkt)

        # Move to the next feature
        feature = layer.GetNextFeature()

    # Save the KML to a file
    #sFilename_kml = 'output.kml'
    kml.save(sFilename_kml)
    return