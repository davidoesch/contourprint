from osgeo import ogr, osr
import numpy as np

# Open the input file geodatabase
input_ds = ogr.Open('/home/menas/Downloads/swissbuilinggdb/swissBUILDINGS3D_3-0_1166-42.gdb', 0)
input_layer_roof = input_ds.GetLayer("Roof_solid")

# Create output geopackage
driver = ogr.GetDriverByName("GPKG")
output_ds = driver.CreateDataSource("roof_points.gpkg")

# Get spatial reference from input layer
spatial_ref = input_layer_roof.GetSpatialRef()

# Create output layer for points
output_layer = output_ds.CreateLayer("roof_points", spatial_ref, ogr.wkbPoint25D)

# Add fields for attributes
field_defn = ogr.FieldDefn("point_id", ogr.OFTInteger)
output_layer.CreateField(field_defn)

# Create feature definition for points
feature_defn = output_layer.GetLayerDefn()

def add_point(x, y, z, point_id):
    point_feature = ogr.Feature(feature_defn)
    point_geom = ogr.Geometry(ogr.wkbPoint25D)
    point_geom.AddPoint(x, y, z)
    point_feature.SetGeometry(point_geom)
    point_feature.SetField("point_id", point_id)
    output_layer.CreateFeature(point_feature)
    return point_id + 1

def handle_geometry(geom, point_id):
    # If geometry has sub-geometries, process them recursively
    if geom.GetGeometryCount() > 0:
        for i in range(geom.GetGeometryCount()):
            subgeom = geom.GetGeometryRef(i)
            point_id = handle_geometry(subgeom, point_id)
    else:
        # Extract points directly from the geometry
        for i in range(geom.GetPointCount()):
            point = geom.GetPoint(i)
            point_id = add_point(point[0], point[1], point[2], point_id)

    return point_id

point_id = 0

boundary_ds = ogr.Open('/home/menas/Downloads/swissbuilinggdb/rectangle.gpkg', 0)
boundary_layer = boundary_ds.GetLayer()
# Create a geometry from the boundary extent
min_x, max_x, min_y, max_y = boundary_layer.GetExtent()
ring = ogr.Geometry(ogr.wkbLinearRing)
ring.AddPoint(min_x, min_y)
ring.AddPoint(max_x, min_y)
ring.AddPoint(max_x, max_y)
ring.AddPoint(min_x, max_y)
ring.AddPoint(min_x, min_y)
boundary_geom = ogr.Geometry(ogr.wkbPolygon)
boundary_geom.AddGeometry(ring)

#list all partial roofs of UUID  which are within the boundary
objectid_list = []
for feature in input_layer_roof:
    geom = feature.GetGeometryRef()
    if boundary_geom.Contains(geom):
        objectid_list.append(feature.GetField("UUID"))

print(len(objectid_list))
counter=0
# Loop through all features wthin boundary and extract information
for objectid in objectid_list:
    print(f"working on roof {counter} of {len(objectid_list)}")
    # Get the feature with the specific EGID in the input layer
    input_layer_roof.SetAttributeFilter(f"UUID= '{objectid}'")
    feature_roof = input_layer_roof.GetNextFeature()

    #check if feature_roof.GetField("EGID") is not None
    if feature_roof.GetField("EGID") is not None:

            geometry = feature_roof.GetGeometryRef()

            # Check if geometry is valid
            if geometry is None:
                continue

            geom_type = geometry.GetGeometryName()
            print(f"Processing Geometry type: {geom_type}")

            # Process the geometry recursively
            point_id = handle_geometry(geometry, point_id)

            print(f"Points processed so far: {point_id}")
    counter=counter+1
print(f"Total points extracted: {point_id}")

# Clean up
input_ds = None
output_ds = None