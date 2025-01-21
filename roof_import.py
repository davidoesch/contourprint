from osgeo import ogr, osr
import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib_scalebar.scalebar import ScaleBar
from owslib.wms import WebMapService
from PIL import Image
from io import BytesIO
import geopandas as gpd
from shapely.geometry import MultiLineString, Polygon
import pandas as pd
import numpy as np
from shapely.ops import unary_union
import math
import subprocess
import contextlib
import io
import warnings

def save_boundary_geometry(geometry, roof_height_max,roof_height_min, uuid, egid, output_path):

    """
    Save boundary polygon to GeoPackage

    Args:
        geometry: shapely Polygon (the boundary_polygon)
        height_type: int
        height: float
        egid: int
        output_path: str, path to output GeoPackage file
    """
    try:
        # Create GeoDataFrame
        gdf = gpd.GeoDataFrame(
            {
                'geometry': [geometry],
                'MAX_HEIGHT': [roof_height_max],
                'MIN_HEIGHT': [roof_height_min],
                'UUID': [uuid],
                'EGID': [egid],
                'area': [geometry.area]
            },
            crs="EPSG:2056"  # Swiss coordinate system
        )
        layer_name = 'height'  # or whatever layer name you prefer

        # Save to GeoPackage
        if os.path.exists(output_path):
            try:
                # Try to read existing layer
                existing_gdf = gpd.read_file(output_path, layer=layer_name)
                # Append the new data
                gdf = pd.concat([existing_gdf, gdf], ignore_index=True)
                # Save with mode='w' to overwrite existing layer
                gdf.to_file(output_path, layer=layer_name, driver="GPKG", mode='w')
            except ValueError:
                # If layer doesn't exist, create new layer
                gdf.to_file(output_path, layer=layer_name, driver="GPKG", mode='a')
        else:
            # If file doesn't exist, create new file
            gdf.to_file(output_path, layer=layer_name, driver="GPKG")

        #print(f"Successfully saved height to {output_path}, layer: {layer_name}")

    except Exception as e:
        print(f"An error occurred while saving height: {str(e)}")

def get_outer_boundary(geom):
    geom_type = geom.GetGeometryType()
    if geom_type != 1016:  # Direct check for TIN type
        raise ValueError(f"Expected TIN geometry, got: {geom_type}")

    # Collect all points from all triangles and create small buffers
    buffers = []
    for i in range(geom.GetGeometryCount()):
        triangle = geom.GetGeometryRef(i)
        ring = triangle.GetGeometryRef(0)

        points = []
        for j in range(ring.GetPointCount()):
            point = ring.GetPoint(j)
            points.append((point[0], point[1]))

        # Create polygon from triangle points and buffer it
        triangle_poly = Polygon(points)
        buffered_triangle = triangle_poly.buffer(0.1)  # 0.1m buffer
        buffers.append(buffered_triangle)

    # Merge all buffers
    merged_buffer = unary_union(buffers)

    return merged_buffer
def convert_gpkg_to_shp(gpkg_path):
    """
    Convert a GeoPackage file to SHP format.

    Args:
        gpkg_path: str, path to the input GeoPackage file
        layer_name: str, name of the layer to convert
        output_path: str, path to the output SHP file
    """
    output_path = gpkg_path.replace('.gpkg', '.shp')
    command = [
        'ogr2ogr',
        '-f', 'ESRI Shapefile',  # Output format
        output_path,  # Output file path
        gpkg_path,  # Input file path

    ]

    try:
        subprocess.run(command, check=True)
        print(f"Successfully exported to {output_path}")
    except subprocess.CalledProcessError as e:
        print(f"Error exporting to SHP: {e}")



# Enable exceptions
ogr.UseExceptions()

# Open the input files
root_url=("/home/menas/Downloads/swissbuilinggdb/")
input_ds = ogr.Open('/home/menas/Downloads/swissbuilinggdb/swissBUILDINGS3D_3-0_1166-42.gdb', 0)
boundary_ds = ogr.Open('/home/menas/Downloads/swissbuilinggdb/rectangle.gpkg', 0)
output_height_geometry = "/home/menas/Downloads/swissbuilinggdb/heightlines.gpkg"
# Check if the output file exists, if yes, delete it
if os.path.exists(output_height_geometry):
    os.remove(output_height_geometry)

# geT layers

boundary_layer = boundary_ds.GetLayer()
input_layer = input_ds.GetLayer("Building_solid")
input_layer_roof = input_ds.GetLayer("Roof_solid")



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
breakpoint()
counter=0
# Loop through all features wthin boundary and extract information
for objectid in objectid_list:
    print(f"working on roof {counter} of {len(objectid_list)}")
    # Get the feature with the specific EGID in the input layer
    input_layer_roof.SetAttributeFilter(f"UUID= '{objectid}'")
    feature_roof = input_layer_roof.GetNextFeature()

    #check if feature_roof.GetField("EGID") is not None
    if feature_roof.GetField("EGID") is not None:

        #get variables ROOF
        roof_max = round(feature_roof.GetField("DACH_MAX"),2)
        roof_min = round(feature_roof.GetField("DACH_MIN"),2)
        egid = feature_roof.GetField("EGID")

        #get variables BUILDING
        input_layer.SetAttributeFilter(f"EGID = '{egid}'")
        feature = input_layer.GetNextFeature()
        gesamthoehe = round(feature.GetField("GESAMTHOEHE"), 2)
        b_dach_max = feature.GetField("DACH_MAX")
        b_dach_min = feature.GetField("DACH_MIN")
        b_gelaendehoehe = feature.GetField("GELAENDEPUNKT")

        #calculate roof height
        # roof_height_max = round(roof_max- b_gelaendehoehe, 2)
        # roof_height_min = round(roof_min- b_gelaendehoehe, 2)

        #calculate roof altitude above ground
        roof_height_max = roof_max
        roof_height_min = roof_min

        # Get the geometry of the feature
        geom_roof = feature_roof.GetGeometryRef()
        outer_ring = get_outer_boundary(geom_roof)

        save_boundary_geometry(geometry=outer_ring , roof_height_max =roof_height_max , roof_height_min=roof_height_min , uuid=objectid, egid=egid, output_path=output_height_geometry)


    counter += 1


# Convert the output GeoPackage to Shapefile
convert_gpkg_to_shp(output_height_geometry)

# Fetch WMS map
bg_layers = ['ch.kantone.cadastralwebmap-farbe', 'ch.swisstopo.swissimage']

for bg_layer in bg_layers:
    print(f"Fetching WMS map for layer: {bg_layer}")
    wms_url = 'https://wms.geo.admin.ch/'
    layer = bg_layer
    wms = WebMapService(wms_url, version='1.3.0')

    # Calculate width and height based on the boundary extent and scale of 1:500
    min_x, max_x, min_y, max_y = boundary_geom.GetEnvelope()
    scale = 500  # 1:500 scale
    # Convert map units to pixels (1 meter = 2 pixels for higher resolution)
    pixels_per_meter = 2
    width = int((max_x - min_x) * pixels_per_meter)
    height = int((max_y - min_y) * pixels_per_meter)

    # Suppress all warnings
    warnings.filterwarnings("ignore")

    with contextlib.redirect_stdout(io.StringIO()): #Suppress warnings
        response = wms.getmap(
            layers=[layer],
            srs='EPSG:2056',
            bbox=(min_x, min_y, max_x, max_y),
            size=(width, height),
            format='image/png',
            transparent=True
        )
    warnings.filterwarnings("default")

    # Create figure with A0 dimensions
    a0_width_inches = 33.1
    a0_height_inches = 23.4
    fig = plt.figure(figsize=(a0_width_inches, a0_height_inches), dpi=300)

    # Calculate the map dimensions in inches at 1:500 scale
    # 1 meter = 0.0394 inches, then divide by scale factor
    map_width_inches = (max_x - min_x) * 0.0394 / scale
    map_height_inches = (max_y - min_y) * 0.0394 / scale

    # Calculate margins to center the map
    left_margin = (a0_width_inches - map_width_inches) / 2
    bottom_margin = (a0_height_inches - map_height_inches) / 2

    # Create axes with the correct size and position
    ax = fig.add_axes([
        left_margin / a0_width_inches,  # left
        bottom_margin / a0_height_inches,  # bottom
        map_width_inches / a0_width_inches,  # width
        map_height_inches / a0_height_inches  # height
    ])

    # Plot the map
    img = Image.open(BytesIO(response.read()))
    ax.imshow(img, extent=[min_x, max_x, min_y, max_y], alpha=0.75)

    # Set limits and aspect ratio
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_y, max_y)
    ax.set_aspect('equal')

    # Format axis labels
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda val, pos: f'{int(round(val))}'))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda val, pos: f'{int(round(val))}'))

    # Add scale bar
    ax.add_artist(ScaleBar(1))

    # Plot geometries from GeoPackage
    gdf = gpd.read_file(output_height_geometry)
    for idx, row in gdf.iterrows():
        polygon = row['geometry']
        max_height = row['MAX_HEIGHT']
        min_height = row['MIN_HEIGHT']

        # Plot polygon
        x, y = polygon.exterior.xy
        ax.plot(x, y, color='blue', linewidth=0.5)
        centroid = polygon.centroid

        # Add height annotations
        for i, point in enumerate(polygon.exterior.coords[:-1]):
            next_point = polygon.exterior.coords[i + 1]
            mid_x = (point[0] + next_point[0]) / 2
            mid_y = (point[1] + next_point[1]) / 2
            angle = math.degrees(math.atan2(next_point[1] - point[1], next_point[0] - point[0]))
            length = np.sqrt((next_point[0] - point[0])**2 + (next_point[1] - point[1])**2)

            if length > 5:
                label_x = point[0] + 0.25 * (next_point[0] - point[0])
                label_y = point[1] + 0.25 * (next_point[1] - point[1])

                ax.text(label_x, label_y, f'{min_height}', color='red', fontsize=3,
                       ha='center', va='center', rotation=angle, rotation_mode='anchor',
                       bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=0.1))

        # Add max height annotation
        ax.text(centroid.x, centroid.y, f'{max_height}', color='green', fontsize=4,
               ha='center', va='center',
               bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=0.1))

    # Save the map
    maptype = bg_layer.split('.')[-1]
    plt.savefig(root_url + "heightinfo_map_" + maptype + ".png", bbox_inches='tight', dpi=300)
    plt.close()

# Clean up
input_ds = None
boundary_ds = None
output_ds = None