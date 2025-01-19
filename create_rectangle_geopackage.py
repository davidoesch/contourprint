import argparse
import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon,MultiPolygon,MultiLineString
from shapely.geometry import LineString, Point
from shapely.geometry import box
from shapely.ops import polygonize
from shapely.ops import unary_union
import numpy as np
import re
import os
import requests
from pyproj import Transformer
import rasterio
from rasterio.merge import merge
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.mask import mask
import subprocess
import sys
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib_scalebar.scalebar import ScaleBar
import math
from osgeo import ogr, osr
from owslib.wms import WebMapService
from PIL import Image
from io import BytesIO
import contextlib
import io
import warnings
import zipfile

"""
create_rectangle_geopackage.py

This script processes geospatial data to generate contour lines, smooth geometries, and extract roof heights.
It also creates maps with height information and exports the results to various formats including PDF, DXF, and PNG.

Usage:
    .venv/bin/python create_rectangle_geopackage.py --bbox "2'604’333/1'197'467 2'604’543/1'197'467 2'604’543/1'197'317 2'604’333/1'197'317" --interval 0.20

Arguments:
    --bbox (str): Bounding box coordinates in the format X/Y, separated by spaces.
                  Default: "2'604’333/1'197'467 2'604’543/1'197'467 2'604’543/1'197'317 2'604’333/1'197'317"
    --interval (float): Interval in meters for contour lines (äquidisztanz). Default is 0.2m.
    --resampling (float): Resampling in meters. Default is 10.0.
    --scale (float): Scale in 1:x in meters. Default is 500.

Functions:
    save_boundary_geometry(geometry, roof_height_max, roof_height_min, uuid, egid, output_path):
        Save boundary polygon to GeoPackage.

    get_outer_boundary(geom):
        Get the outer boundary of a TIN geometry.

    convert_gpkg_to_shp(gpkg_path):
        Convert a GeoPackage file to SHP format.

    parse_coordinates(coord_str):
        Parse a coordinate string into a tuple of integers.

    create_rectangle_geopackages(coord_string):
        Create GeoPackages with a rectangle polygon and a 10-meter buffer.

    download_swissalti3d(buffered_rectangle_gpkg):
        Download the latest swissALTI3D data intersecting the buffered rectangle.

    download_swissbuilding3d(buffered_rectangle_gpkg):
        Download the latest Swissbuilding data intersecting the buffered rectangle.

    clip_and_resample_dem(buffered_rectangle_gpkg, input_tif, output_tif="swissalti_clipped.tif", target_resolution=None):
        Clip and optionally resample the input DEM TIFF file.

    generate_contours(interval):
        Generate contour lines from a raster using gdal_contour.

    convert_close_contours_to_polygons(gdf, tolerance=1e-6):
        Convert contour lines to polygons if their start and end points are close.

    process_contours(interval):
        Process contours by generating and converting them to polygons.

    smooth_geometry(interval, input_file, output_file):
        Smooth the geometry of the lines in the specified input file using QGIS processing.

    merge_smoothed():
        Merge smoothed contour polygons and lines into a single GeoPackage.

    export_to_image_pdf(interval):
        Export the contours and rectangle to PDF, PNG, and DXF formats.

    check_qgis_installed():
        Check if QGIS is installed and accessible.

    process_roof():
        Process roof data to extract height information and save it to a GeoPackage.

    create_map(maptype, scale, interval):
        Create a map with the specified type, scale, and interval.


"""

# example call: .venv/bin/python create_rectangle_geopackage.py --bbox "2'604’333/1'197'467 2'604’543/1'197'467 2'604’543/1'197'317 2'604’333/1'197'317" --interval 0.20

global resample
global output_folder
destionation_name = "results"
output_folder = os.path.join(os.getcwd(), destionation_name )
BUFFER_M = 10
SWISSALTI3D_COLLECTION = "ch.swisstopo.swissalti3d"
SWISSBUILDING3D_COLLECTION = "ch.swisstopo.swissbuildings3d_3_0"
STAC_API_URL = "https://data.geo.admin.ch/api/stac/v0.9"



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

def parse_coordinates(coord_str):
    """
    Parses a coordinate string in the format "2'604’333/1'197'467" into a tuple of integers.
    """
    x_str, y_str = re.sub(r"[^\d/]", "", coord_str).split("/")
    return int(x_str), int(y_str)

def create_rectangle_geopackages(coord_string):

    """
    Creates two GeoPackages: one with a rectangle polygon and one with a 10-meter buffer.

    """
    os.makedirs(output_folder, exist_ok=True)
    rectangle_output=os.path.join(output_folder,"rectangle.gpkg")

    buffer_output=os.path.join(output_folder,"rectangle_buffer.gpkg")


    coord_strings = coord_string.split()
    if len(coord_strings) != 4:
        raise ValueError("Please provide exactly 4 coordinate pairs for the rectangle.")
    coordinates = [parse_coordinates(coord) for coord in coord_strings]

    # Create rectangle and buffer geometries
    rectangle = Polygon(coordinates)
    gdf = gpd.GeoDataFrame([{'geometry': rectangle}], crs="EPSG:2056")
    buffer_gdf = gpd.GeoDataFrame([{'geometry': rectangle.buffer(BUFFER_M)}], crs="EPSG:2056")

    # Save to GeoPackages
    gdf.to_file(rectangle_output, layer="rectangle_layer", driver="GPKG")
    buffer_gdf.to_file(buffer_output, layer="buffer_layer", driver="GPKG")
    print(f"Rectangle GeoPackage saved to {rectangle_output}")
    print(f"Buffered GeoPackage saved to {buffer_output}")

    return buffer_output  # Return the path to the buffered GeoPackage


def download_swissalti3d(buffered_rectangle_gpkg):
    """
    Downloads the latest swissALTI3D data intersecting the buffered rectangle.
    The bounding box is transformed from EPSG:2056 to EPSG:4326.
    """
    # Initialize transformer from EPSG:2056 to EPSG:4326
    transformer = Transformer.from_crs("EPSG:2056", "EPSG:4326", always_xy=True)

    # Read the buffered rectangle geometry from GeoPackage
    buffered_gdf = gpd.read_file(buffered_rectangle_gpkg)

    # Calculate bounding box of buffered area
    minx, miny, maxx, maxy = buffered_gdf.total_bounds

    # Transform the bounding box to WGS84
    min_lon, min_lat = transformer.transform(minx, miny)
    max_lon, max_lat = transformer.transform(maxx, maxy)
    bbox = [min_lon, min_lat, max_lon, max_lat]

    # Query STAC API for items intersecting the bounding box
    items_url = f"{STAC_API_URL}/collections/{SWISSALTI3D_COLLECTION}/items?bbox={','.join(map(str, bbox))}&limit=100&sortby=-datetime"
    print(f"Trying to connect to STAC API: ")
    try:
        response = requests.get(items_url,timeout=15)
        items = response.json().get("features", [])
    except requests.exceptions.RequestException as e:
        print(f"Failed to connect to STAC API: {e}")
        return

    if not items:
        print("No matching data found for the specified area.")
        return

    # Ensure output folder exists
    os.makedirs(output_folder, exist_ok=True)

    # Download assets for each matching item
    tiff_files = []  # List to store paths of downloaded TIFF files
    for item in items:
        item_id = item["id"]
        assets = item.get("assets", {})

        # Filter assets based on the criteria
        for asset_key, asset in assets.items():
            if "_0.5_" in asset_key and asset_key.endswith(".tif"):
                asset_url = asset["href"]
                asset_filename = os.path.join(output_folder, f"{item_id}_{asset_key}")

                print(f"Downloading {asset_url} to {asset_filename}")
                try:
                    asset_data = requests.get(asset_url,timeout=15)
                    asset_data.raise_for_status()  # Raise an error for bad responses

                    with open(asset_filename, "wb") as f:
                        for chunk in asset_data.iter_content(chunk_size=8192):
                            f.write(chunk)

                    # Add the TIFF file to the list for merging
                    tiff_files.append(asset_filename)

                except requests.exceptions.RequestException as e:
                    print(f"Failed to download {asset_url}: {e}")

    print(f"Downloaded data is saved in '{output_folder}'")

    # Define the output filename
    output_filename = os.path.join(output_folder, "swissalti.tif")

    # Merge TIFF files if more than one is found
    if len(tiff_files) > 1:

        # Create a list to store the opened raster files
        src_files_to_mosaic = []

        # Open each TIFF file
        for fp in tiff_files:
            src = rasterio.open(fp)
            src_files_to_mosaic.append(src)

        # Merge the rasters
        mosaic, out_trans = merge(src_files_to_mosaic)

        # Copy the metadata from one of the source files
        out_meta = src_files_to_mosaic[0].meta.copy()

        # Update the metadata
        out_meta.update({
            "driver": "GTiff",
            "height": mosaic.shape[1],
            "width": mosaic.shape[2],
            "transform": out_trans
        })

        # Write the mosaic to a new TIFF file
        with rasterio.open(output_filename, "w", **out_meta) as dest:
            dest.write(mosaic)

        print(f"Merged TIFF files into '{output_filename}'")
    elif tiff_files:
        # If only one TIFF file found, just copy it to the new filename
        single_file = tiff_files[0]
        os.rename(single_file, output_filename)
        print(f"Single TIFF file saved as '{output_filename}'")
    else:
        print("No TIFF files downloaded.")

    # Clip the swissalti.tif with the buffered rectangle
    clip_and_resample_dem(buffered_rectangle_gpkg, output_filename,target_resolution=resample)

def download_swissbuilding3d(buffered_rectangle_gpkg):
    """
    Downloads the latest swSwissbuilding data intersecting the buffered rectangle.
    The bounding box is transformed from EPSG:2056 to EPSG:4326.
    """
    # Initialize transformer from EPSG:2056 to EPSG:4326
    transformer = Transformer.from_crs("EPSG:2056", "EPSG:4326", always_xy=True)

    # Read the buffered rectangle geometry from GeoPackage
    buffered_gdf = gpd.read_file(buffered_rectangle_gpkg)

    # Calculate bounding box of buffered area
    minx, miny, maxx, maxy = buffered_gdf.total_bounds

    # Transform the bounding box to WGS84
    min_lon, min_lat = transformer.transform(minx, miny)
    max_lon, max_lat = transformer.transform(maxx, maxy)
    bbox = [min_lon, min_lat, max_lon, max_lat]

    # Query STAC API for items intersecting the bounding box
    items_url = f"{STAC_API_URL}/collections/{SWISSBUILDING3D_COLLECTION}/items?bbox={','.join(map(str, bbox))}&limit=100&sortby=-datetime"
    print(f"Trying to connect to STAC API: ")
    try:
        response = requests.get(items_url,timeout=15)
        items = response.json().get("features", [])
    except requests.exceptions.RequestException as e:
        print(f"Failed to connect to STAC API: {e}")
        return

    if not items:
        print("No matching data found for the specified area.")
        return

    # Ensure output folder exists
    os.makedirs(output_folder, exist_ok=True)

    # Download assets for each matching item
    GDB_files = []  # List to store paths of downloaded TIFF files
    for item in items:
        item_id = item["id"]
        assets = item.get("assets", {})

        # Filter assets based on the criteria
        for asset_key, asset in assets.items():

            if ".gdb." in asset_key and "_2021_" in asset_key and asset_key.endswith(".zip"):
                print("only downloading 2021 data")
                asset_url = asset["href"]
                asset_filename = os.path.join(output_folder, f"{item_id}_{asset_key}")

                print(f"Downloading {asset_url} to {asset_filename}")
                try:
                    asset_data = requests.get(asset_url,timeout=15)
                    asset_data.raise_for_status()  # Raise an error for bad responses

                    with open(asset_filename, "wb") as f:
                        for chunk in asset_data.iter_content(chunk_size=8192):
                            f.write(chunk)

                    # Add the TIFF file to the list for merging
                    GDB_files.append(asset_filename)

                except requests.exceptions.RequestException as e:
                    print(f"Failed to download {asset_url}: {e}")

    print(f"Downloaded data is saved in '{output_folder}'")

    # Define the output filename
    output_filename = os.path.join(output_folder, "swissbuilding3D.gdb")
    # Merge TIFF files if more than one is found
    if len(GDB_files) > 1:


        print(f"Multiple GDB files, canot be merged, files stoired here '{output_filename}': adapt the code: give back the list of files, then in roof_process use the list to loop trough the files")
    elif GDB_files:
        # If only one TIFF file found, just copy it to the new filename
        single_file = GDB_files[0]
        # unzip singel file to output folder with output_filename
        with zipfile.ZipFile(single_file, 'r') as zip_ref:
            temp_extract_path = os.path.join(output_folder, "temp_extract")
            zip_ref.extractall(temp_extract_path)
            # Move the extracted folder to the desired output path
            extracted_folder_name = os.listdir(temp_extract_path)[0]
            extracted_folder_path = os.path.join(temp_extract_path, extracted_folder_name)
            if os.path.exists(output_filename):
                if os.path.isdir(output_filename):
                    import shutil
                    shutil.rmtree(output_filename)
                else:
                    os.remove(output_filename)
            os.rename(extracted_folder_path, output_filename)
            # Clean up the temporary extraction path
            os.rmdir(temp_extract_path)
        #os.rename(single_file, output_filename)

        print(f"Single GDB file saved as '{output_filename}'")
    else:
        print("No GDB files downloaded.")



def clip_and_resample_dem(buffered_rectangle_gpkg, input_tif, output_tif="swissalti_clipped.tif", target_resolution=None):
    """
    Clips the input TIFF using the buffered rectangle geometry from a GeoPackage,
    optionally resamples to a target resolution with an additional buffer to mitigate edge effects,
    and saves the result.

    :param buffered_rectangle_gpkg: Path to the GeoPackage containing the buffered rectangle
    :param input_tif: Path to the input DEM TIFF file
    :param output_tif: Name of the output TIFF file (default: "swissalti_clipped.tif")
    :param target_resolution: Target resolution in meters (e.g., 5 for 5m resolution). If None, no resampling is performed.
    """
    # Read the buffered rectangle from the GeoPackage
    buffered_gdf = gpd.read_file(buffered_rectangle_gpkg)

    # Assuming there's only one geometry in the buffered GeoDataFrame
    buffered_geom = buffered_gdf.geometry.values[0]

    # Create a GeoJSON-like dictionary for the geometry
    geojson_geom = {
        "type": "Polygon",
        "coordinates": [list(buffered_geom.exterior.coords)]
    }

    # Open the input TIFF file
    with rasterio.open(input_tif) as src:
        if target_resolution is not None:
            # Create an additional buffer of 100m for resampling
            additional_buffer = 100  # meters
            extended_bounds = box(*buffered_geom.bounds).buffer(additional_buffer)
            extended_geojson = {
                "type": "Polygon",
                "coordinates": [list(extended_bounds.exterior.coords)]
            }

            # Clip the raster with the extended mask
            out_image, out_transform = mask(src, [extended_geojson], crop=True)
        else:
            # If no resampling, use the original buffered geometry
            out_image, out_transform = mask(src, [geojson_geom], crop=True)

        # Check the shape of out_image and squeeze if necessary
        if out_image.ndim == 3:  # Shape (1, height, width)
            out_image = out_image.squeeze()  # Remove extra dimension

        # Get the metadata and update it
        out_meta = src.meta.copy()
        out_meta.update({
            "driver": "GTiff",
            "height": out_image.shape[0],
            "width": out_image.shape[1],
            "transform": out_transform,
        })

        # Resample if target_resolution is specified
        if target_resolution is not None:
            # Calculate scaling factor
            current_resolution = src.res[0]  # Assuming square pixels
            scale_factor = current_resolution / target_resolution

            # Calculate new dimensions
            new_height = int(out_image.shape[0] * scale_factor)
            new_width = int(out_image.shape[1] * scale_factor)

            # Update transform for the new resolution
            new_transform = rasterio.transform.from_bounds(
                *rasterio.transform.array_bounds(out_image.shape[0], out_image.shape[1], out_transform),
                new_width,
                new_height
            )

            # Perform resampling
            resampled_image = np.zeros((1, new_height, new_width), dtype=out_image.dtype)
            reproject(
                source=out_image[np.newaxis, :, :],
                destination=resampled_image,
                src_transform=out_transform,
                src_crs=src.crs,
                dst_transform=new_transform,
                dst_crs=src.crs,
                resampling=Resampling.bilinear
            )

            # Update image and metadata
            out_image = resampled_image.squeeze()
            out_meta.update({
                "height": new_height,
                "width": new_width,
                "transform": new_transform,
            })

            # Write the resampled image to a temporary file
            temp_tif = "temp_resampled.tif"
            with rasterio.open(temp_tif, "w", **out_meta) as temp_dst:
                temp_dst.write(out_image, 1)

            # Clip the resampled image back to the original buffered rectangle
            with rasterio.open(temp_tif) as temp_src:
                out_image, out_transform = mask(temp_src, [geojson_geom], crop=True)

            out_image = out_image.squeeze()

            # Update metadata again after final clipping
            out_meta.update({
                "height": out_image.shape[0],
                "width": out_image.shape[1],
                "transform": out_transform,
            })

            # Remove the temporary file
            os.remove(temp_tif)

        # Save the clipped and optionally resampled raster to a new file
        output_tif_path = os.path.join(os.path.dirname(input_tif), output_tif)
        with rasterio.open(output_tif_path, "w", **out_meta) as dst:
            dst.write(out_image, 1)  # Write the single band

    print(f"Clipped{'and resampled ' if target_resolution else ' '}TIFF file saved as '{output_tif}'")


def generate_contours(interval):
    """
    Generates contour lines from a raster using gdal_contour.

    Args:
        interval (float): Contour interval (default is 0.5).
        band (int): Band number to use for contouring (default is 1).
        attribute (str): Attribute name to assign elevation values (default is "ELEV").
    """
    band=1
    attribute="ELEV"
    input_tif = os.path.join(output_folder, "swissalti_clipped.tif")
    output_gpkg = os.path.join(output_folder, "contour.gpkg")

    command = [
        "gdal_contour",
        f"-b", str(band),          # Band number
        f"-a", attribute,          # Attribute name for elevation
        f"-i", str(interval),      # Interval
        f"-f", "GPKG",             # Output format
        input_tif,                # Input raster file
        output_gpkg                # Output GeoPackage
    ]

    try:
        subprocess.run(command)#, check=True)
        print(f"Contour lines saved to '{output_gpkg}'")
    except subprocess.CalledProcessError as e:
        print(f"Error generating contours: {e}")

def convert_close_contours_to_polygons(gdf, tolerance=1e-6):
    """
    Converts contour lines to polygons if their start and end points are close.

    Args:
        gdf (GeoDataFrame): Input GeoDataFrame containing contour lines.
        tolerance (float): Maximum distance between start and end points to consider as closed.

    Returns:
        GeoDataFrame: A new GeoDataFrame with polygons where applicable.
    """
    def line_to_polygon(geom):
        if isinstance(geom, LineString):
            if geom.is_closed or Point(geom.coords[0]).distance(Point(geom.coords[-1])) <= tolerance:
                return Polygon(geom)
        elif isinstance(geom, Polygon):
            return geom
        return geom  # Return original geometry if not LineString or Polygon

    # Apply the conversion
    gdf['geometry'] = gdf['geometry'].apply(line_to_polygon)

    # Filter to keep only Polygon and MultiPolygon geometries
    return gdf[gdf['geometry'].type.isin(['Polygon', 'MultiPolygon'])]

def process_contours(interval):
    generate_contours(interval)

    # Read the generated contours
    contours_gpkg = os.path.join(output_folder, "contour.gpkg")
    gdf = gpd.read_file(contours_gpkg)

    # Convert close contours to polygons
    polygons_gdf = convert_close_contours_to_polygons(gdf)

    # Save the polygons
    polygons_output = os.path.join(output_folder, "contour_polygons.gpkg")
    polygons_gdf.to_file(polygons_output, driver="GPKG")
    print(f"Contour polygons saved to '{polygons_output}'")

    # Keep all geometries that are not Polygons or MultiPolygons
    lines_gdf = gdf[~gdf.geometry.apply(lambda geom: isinstance(geom, (Polygon, MultiPolygon)))]

    # Save the remaining geometries (mostly lines)
    lines_output = os.path.join(output_folder, "contour_lines.gpkg")
    lines_gdf.to_file(lines_output, driver="GPKG")
    print(f"Contour lines saved to '{lines_output}'")

def smooth_geometry(interval,input_file,output_file):
    """
    Smooths the geometry of the lines in the specified input file using QGIS processing.

    : interval smoothing minimum distance.  Normally 0.25m
    """

    offset=0.25
    # Define the QGIS process command
    command = [
        "/usr/bin/qgis_process",  # Ensure qgis_process is in your PATH
        "run", "native:smoothgeometry",
        "--distance_units=meters",
        "--area_units=m2",
        "--ellipsoid=EPSG:7004",
        f"--INPUT={input_file}",
        "--ITERATIONS=10",
        f"--OFFSET={offset}",
        "--MAX_ANGLE=180",
        f"--OUTPUT={output_file}"
    ]

    try:
        # Run the command
        result = subprocess.run(command)#, check=True, capture_output=True, text=True)
        print(f"Smoothing {output_file} completed successfully.")
        #print(result.stdout)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred in {output_file} : {e.stderr}")

def merge_smoothed():
    # Read the polygon and line files
    polygons_gdf = gpd.read_file(os.path.join(output_folder, "contour_polygons_smoothed.gpkg"))
    lines_gdf = gpd.read_file(os.path.join(output_folder, "contour_lines_smoothed.gpkg"))

    # Convert polygons to lines
    def polygon_to_line(geometry):
        if geometry.geom_type == 'Polygon':
            return LineString(geometry.exterior.coords)
        elif geometry.geom_type == 'MultiPolygon':
            return MultiLineString([LineString(poly.exterior.coords) for poly in geometry])
        else:
            return geometry

    polygons_as_lines = polygons_gdf.copy()
    polygons_as_lines['geometry'] = polygons_as_lines['geometry'].apply(polygon_to_line)

    # Merge the converted polygons with the existing lines
    merged_gdf = gpd.GeoDataFrame(pd.concat([polygons_as_lines, lines_gdf], ignore_index=True))

    # Save the merged result
    output_file = os.path.join(output_folder, "contour_smoothed.gpkg")
    merged_gdf.to_file(output_file, driver="GPKG")
    print(f"Merged contours saved to '{output_file}'")

def export_to_image_pdf(interval):
    # Step 1: Read the Contour GeoPackage
    gpkg_path = os.path.join(output_folder, "contour_smoothed.gpkg")
    layer_name = 'contour_smoothed'
    output_path = os.path.join(output_folder, "isolinien")
    gdf_contours = gpd.read_file(gpkg_path, layer=layer_name)

    # Ensure the GeoDataFrame is in the correct CRS (EPSG:2056)
    if gdf_contours.crs != 'EPSG:2056':
        gdf_contours = gdf_contours.to_crs(epsg=2056)

    # Step 2: Read the Rectangle GeoPackage
    rectangle_path = os.path.join(output_folder, "rectangle.gpkg")
    gdf_rectangle = gpd.read_file(rectangle_path)

    # Ensure the rectangle GeoDataFrame is in the correct CRS (EPSG:2056)
    if gdf_rectangle.crs != 'EPSG:2056':
        gdf_rectangle = gdf_rectangle.to_crs(epsg=2056)
    #function to set line widths based on elevation
    def get_line_width(elevation):
        if elevation % 10 == 0:
            return 3.5  # Thicker line
        elif elevation % 1 == 0:
            return 2.5  # Medium thickness
        else:
            return 1.0  # Thinner line
    # Apply the function to create a new 'line_width' column
    gdf_contours['line_width'] = gdf_contours['ELEV'].apply(get_line_width)

    # Step 4: Plot the Contours
    min_x, min_y, max_x, max_y = gdf_contours.total_bounds
    aspect_ratio = (max_y - min_y) / (max_x - min_x)

    # Set figure size for A0
    a0_width_in = 841 / 25.4  # Convert mm to inches
    a0_height_in = 1189 / 25.4  # Convert mm to inches
    fig, ax = plt.subplots(figsize=(a0_width_in, a0_height_in), dpi=300)

    # Plot each contour line
    for _, row in gdf_contours.iterrows():
        if row['geometry'] is not None:
            x, y = row['geometry'].xy
            ax.plot(x, y, linewidth=row['line_width'], color='black')

            # Label the thickest lines (10m increment lines)
            if row['line_width'] == 3.5:
                line = LineString(row['geometry'])
                if line.length < 100:
                    point = line.interpolate(line.length / 2)
                    x_label, y_label = point.x, point.y
                    angle = 0
                    ax.text(x_label, y_label, f"{int(row['ELEV'])} m", fontsize=16, ha='center', va='center', color='blue', rotation=angle, rotation_mode='anchor')
                else:
                    num_labels = int(line.length // 60)
                    for i in range(num_labels):
                        distance = i * 60
                        point = line.interpolate(distance)
                        x_label, y_label = point.x, point.y
                        if distance + 1 < line.length:
                            next_point = line.interpolate(distance + 1)
                            dx = next_point.x - x_label
                            dy = next_point.y - y_label
                            angle = np.degrees(np.arctan2(dy, dx))
                        else:
                            angle = 0
                        ax.text(x_label, y_label, f"{int(row['ELEV'])} m", fontsize=16, ha='center', va='center', color='blue', rotation=angle, rotation_mode='anchor')

    # Plot the rectangle overlay
    for _, rect_row in gdf_rectangle.iterrows():
        if rect_row['geometry'] is not None:
            x_rect, y_rect = rect_row['geometry'].exterior.xy
            ax.plot(x_rect, y_rect, color='red', linewidth=2)

    # Step 5: Set Axis Limits Based on GeoPackage Extent
    ax.set_xlim(min_x, max_x)
    ax.set_ylim(min_y, max_y)

    # Keep the aspect ratio of the plot
    ax.set_aspect('equal', adjustable='box')

    # Adjust ticks to display scale
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda val, pos: f'{int(round(val))}'))
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda val, pos: f'{int(round(val))}'))

    # Step 6: Add a Simple Scale Bar
    ax.add_artist(ScaleBar(1))

    # Add a white rectangle behind the label
    rect = patches.Rectangle((0.001, 0.001), 0.05, 0.01, transform=ax.transAxes, color='white', zorder=2)
    ax.add_patch(rect)

    # Add the 'Äquidistanz' and 'DEM' label
    display_value = 0.5 if resample is None else resample
    ax.text(0.001, 0.001, 'Äquidistanz: ' + str(interval)+'m'+' DEM: '+str(display_value)+'m ' , ha='left', va='bottom', transform=ax.transAxes,
            fontsize=14, color='black', zorder=3)

    # Step 7: Save as Image or PDF
    plt.savefig(output_path + '.pdf', bbox_inches='tight', dpi=300)
    plt.savefig(output_path + '.png', bbox_inches='tight', dpi=300)

    # Step 8: Export the contours and rectahngel to DXF using ogr2ogr
    dxf_output_path = output_path + '.dxf'
    command = [
        'ogr2ogr',
        '-f', 'DXF',  # Output format
        dxf_output_path,  # Output file path
        gpkg_path,  # Input file path
        layer_name  # Layer name to export
    ]

    try:
        subprocess.run(command, check=True)
        print(f"Successfully exported to {dxf_output_path}")
    except subprocess.CalledProcessError as e:
        print(f"Error exporting to DXF: {e}")

    rect_output_path = os.path.join(output_folder,'rectangle.dxf')
    command = [
        'ogr2ogr',
        '-f', 'DXF',  # Output format
        rect_output_path,  # Output file path
        rectangle_path,  # Input file path
        "rectangle_layer"  # Layer name to export
    ]

    try:
        subprocess.run(command, check=True)
        print(f"Successfully exported to {rect_output_path}")
    except subprocess.CalledProcessError as e:
        print(f"Error exporting to DXF: {e}")


    plt.close()

def check_qgis_installed():
    try:
        subprocess.run(["/usr/bin/qgis_process", "--version"])#, check=True, capture_output=True, text=True)
        print("QGIS is installed and accessible.")
    except FileNotFoundError:
        print("Warning: QGIS is not installed or not in PATH. This script requires QGIS to run properly.")
        sys.exit(1)  # Exit if QGIS is not found

def process_roof():
    input_ds = ogr.Open(os.path.join(output_folder, "swissbuilding3D.gdb"), 0)
    boundary_ds = ogr.Open(os.path.join(output_folder, "rectangle.gpkg"), 0)
    output_height_geometry = os.path.join(output_folder, "heightlines.gpkg")
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

            roof_height_max = round(roof_max- b_gelaendehoehe, 2)
            roof_height_min = round(roof_min- b_gelaendehoehe, 2)

            # Get the geometry of the feature
            geom_roof = feature_roof.GetGeometryRef()
            outer_ring = get_outer_boundary(geom_roof)

            save_boundary_geometry(geometry=outer_ring , roof_height_max =roof_height_max , roof_height_min=roof_height_min , uuid=objectid, egid=egid, output_path=output_height_geometry)


        counter += 1


    # Convert the output GeoPackage to Shapefile
    convert_gpkg_to_shp(output_height_geometry)

def create_map(maptype,scale,interval):
    # Fetch WMS map
    bg_layers = ['ch.kantone.cadastralwebmap-farbe', 'ch.swisstopo.swissimage']

    for bg_layer in bg_layers:
        print(f"Fetching WMS map for layer: {bg_layer}")
        wms_url = 'https://wms.geo.admin.ch/'
        layer = bg_layer
        wms = WebMapService(wms_url, version='1.3.0')

        # Read the rectangle GeoPackage to get the boundary geometry
        boundary_gdf = gpd.read_file(os.path.join(output_folder, "rectangle.gpkg"))
        boundary_geom = boundary_gdf.geometry.values[0]
        min_x, min_y, max_x, max_y = boundary_geom.bounds

        width = int((max_x - min_x) / (scale/1000))  # 0.5 is 1:500 scale
        height = int((max_y - min_y) / (scale/1000))  # 1:500 scale
        # Suppress all warnings
        warnings.filterwarnings("ignore")

        with contextlib.redirect_stdout(io.StringIO()): #Suppress wanrings
            response = wms.getmap(
                layers=[layer],
                srs='EPSG:2056',
                bbox=(min_x, min_y, max_x, max_y),
                size=(width * 2, height * 2),  # Increase resolution by doubling the size
                format='image/png',
                transparent=True
            )
        # Re-enable warnings after the operation if needed
        warnings.filterwarnings("default")

        # Plot the map and points
        img = Image.open(BytesIO(response.read()))
        plt.figure(figsize=(33.1, 23.4), dpi=300) # A0 size in inches (landscape)
        plt.imshow(img, extent=[min_x, max_x, min_y, max_y], alpha=0.75)


        # Step 5: Set Axis Limits Based on GeoPackage Extent
        plt.xlim(min_x, max_x)
        plt.ylim(min_y, max_y)

        # Keep the aspect ratio of the plot
        plt.gca().set_aspect('equal', adjustable='box')

        # Adjust ticks to display scale
        plt.gca().xaxis.set_major_formatter(plt.FuncFormatter(lambda val, pos: f'{int(round(val))}'))
        plt.gca().yaxis.set_major_formatter(plt.FuncFormatter(lambda val, pos: f'{int(round(val))}'))

        # Add a Simple Scale Bar
        plt.gca().add_artist(ScaleBar(1))

        # open output_height_geometry and plot the polygon. add the height
        if maptype == "heightinfo":
            print("plotting height info")
            gdf = gpd.read_file(os.path.join(output_folder, "heightlines.gpkg"))
            for idx, row in gdf.iterrows():
                polygon = row['geometry']
                max_height = row['MAX_HEIGHT']
                min_height = row['MIN_HEIGHT']

                # Plot the polygon
                x, y = polygon.exterior.xy
                plt.plot(x, y, color='blue', linewidth=0.5)
                centroid = polygon.centroid

                # Annotate the max height on the line of the polygon
                for i, point in enumerate(polygon.exterior.coords[:-1]):
                    next_point = polygon.exterior.coords[i + 1]
                    mid_x = (point[0] + next_point[0]) / 2
                    mid_y = (point[1] + next_point[1]) / 2
                    # Calculate angle of the line segment
                    angle = math.degrees(math.atan2(next_point[1] - point[1], next_point[0] - point[0]))
                    # Calculate length of the line segment
                    length = np.sqrt((next_point[0] - point[0])**2 + (next_point[1] - point[1])**2)

                    # Only add label if the segment is long enough (adjust threshold as needed)
                    if length > 5:  # You may need to adjust this threshold
                        # Calculate position 25% along the line from the start point
                        label_x = point[0] + 0.25 * (next_point[0] - point[0])
                        label_y = point[1] + 0.25 * (next_point[1] - point[1])

                        # Rotate text to align with the line
                        plt.text(label_x, label_y, f'{min_height}', color='red', fontsize=3, ha='center', va='center',
                                rotation=angle, rotation_mode='anchor',
                                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=0.1))



                # Annotate the max height
                plt.text(centroid.x, centroid.y, f'{max_height}', color='green', fontsize=4, ha='center', va='center',
                        bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=0.1))
        if maptype == "isolinien":

            gdf_contours = gpd.read_file(os.path.join(output_folder, "contour_smoothed.gpkg"))



            # Step 2: Read the Rectangle GeoPackage
            rectangle_path = os.path.join(output_folder, "rectangle.gpkg")
            gdf_rectangle = gpd.read_file(rectangle_path)

            # Ensure the rectangle GeoDataFrame is in the correct CRS (EPSG:2056)
            if gdf_rectangle.crs != 'EPSG:2056':
                gdf_rectangle = gdf_rectangle.to_crs(epsg=2056)
            #function to set line widths based on elevation
            def get_line_width(elevation):
                if elevation % 10 == 0:
                    return 3.5  # Thicker line
                elif elevation % 1 == 0:
                    return 2.5  # Medium thickness
                else:
                    return 1.0  # Thinner line
            # Apply the function to create a new 'line_width' column
            gdf_contours['line_width'] = gdf_contours['ELEV'].apply(get_line_width)

            # Plot each contour line
            for _, row in gdf_contours.iterrows():
                if row['geometry'] is not None:
                    x, y = row['geometry'].xy
                    plt.plot(x, y, linewidth=row['line_width'] + 2, color='white')  # White border
                    plt.plot(x, y, linewidth=row['line_width'], color='black')  # Black line

                    # Label every second line
                    if row['line_width'] == 3.5 or row['line_width'] == 2.5:
                        line = LineString(row['geometry'])
                        if line.length < 100:
                            point = line.interpolate(line.length / 2)
                            x_label, y_label = point.x, point.y
                            angle = 0
                            plt.text(x_label, y_label, f"{int(row['ELEV'])} m", fontsize=16, ha='center', va='center', color='blue', rotation=angle, rotation_mode='anchor')
                        else:
                            num_labels = int(line.length // 60)
                            for i in range(num_labels):
                                distance = i * 60
                                point = line.interpolate(distance)
                                x_label, y_label = point.x, point.y
                                if distance + 1 < line.length:
                                    next_point = line.interpolate(distance + 1)
                                    dx = next_point.x - x_label
                                    dy = next_point.y - y_label
                                    angle = np.degrees(np.arctan2(dy, dx))
                                else:
                                    angle = 0
                                plt.text(x_label, y_label, f"{int(row['ELEV'])} m", fontsize=16, ha='center', va='center', color='blue', rotation=angle, rotation_mode='anchor')

            # Plot the rectangle overlay
            for _, rect_row in gdf_rectangle.iterrows():
                if rect_row['geometry'] is not None:
                    x_rect, y_rect = rect_row['geometry'].exterior.xy
                    plt.plot(x_rect, y_rect, color='red', linewidth=2)

            # Add a white rectangle behind the label at the bottom right
            rect = patches.Rectangle((0.85, 0.01), 1.50, 0.03, transform=plt.gca().transAxes, color='white', zorder=2)
            plt.gca().add_patch(rect)

            # Add the 'Äquidistanz' and 'DEM' label at the bottom right
            display_value = 0.5 if resample is None else resample
            plt.text(0.975, 0.015, 'Äquidistanz: ' + str(interval) + 'm' + ' DEM: ' + str(display_value) + 'm',
                     ha='right', va='bottom', transform=plt.gca().transAxes, fontsize=14, color='black', zorder=3)

        # Save the map
        mapname = bg_layer.split('.')[-1]
        plt.savefig(os.path.join(output_folder,maptype+"_map_"+mapname+".png"),bbox_inches='tight', dpi=300)
        plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create Contour Lines for a specified bounding box and interval in EPSG:2056.")

    # BBox argument with flag (optional, defaulting to a standard bounding box)
    parser.add_argument(
        "--bbox",
        type=str,
        default="2'604’333/1'197'467 2'604’543/1'197'467 2'604’543/1'197'317 2'604’333/1'197'317",
        help="Bounding box coordinates in the format X/Y, separated by spaces (default: '2'604’333/1'197'467 2'604’543/1'197'467 2'604’543/1'197'317 2'604’333/1'197'317')."
    )

    # Interval argument with flag (optional, defaulting to 10 meters)
    parser.add_argument(
        "--interval",
        type=float,
        default=.20,
        help="Interval in meters for contour lines (äquidisztanz). Default is 0.2m."
    )

    # Interval argument with flag (optional, defaulting to 10 meters)
    parser.add_argument(
        "--resampling",
        type=float,
        default=10.0,
        help="Resampling in meters . Default is None."
    )

    # Map Scale default 1:500
    parser.add_argument(
        "--scale",
        type=float,
        default=500,
        help="Scale in 1:x in meters . Default is 500."
    )

    args = parser.parse_args()

    print(" -------------------")
    print()
    print("IMPORTANT: DO NOT Run from activated .venv, not from within VSCODE")
    print()
    print(" -------------------")

    # Set the resampling value
    scale=args.scale

    # Set the resampling value
    resample=args.resampling
    output_folder = os.path.join(os.getcwd(), destionation_name ,str(args.resampling))

    # Check for QGIS installation
    check_qgis_installed()


    # Create the contour lines based on input or default arguments
    buffered_gpkg = create_rectangle_geopackages(args.bbox)

    #Download Data
    download_swissalti3d(buffered_gpkg)

    #generate contours
    #generate_contours(interval=args.interval)
    process_contours(interval=args.interval)


    #smooth contours using QGIS
    smooth_geometry(interval=args.interval,
                    input_file=os.path.join(output_folder, "contour_lines.gpkg|layername=contour_lines"),
                    output_file=os.path.join(output_folder, "contour_lines_smoothed.gpkg"))
    smooth_geometry(interval=args.interval,
                    input_file=os.path.join(output_folder, "contour_polygons.gpkg|layername=contour_polygons"),
                    output_file=os.path.join(output_folder, "contour_polygons_smoothed.gpkg"))
    merge_smoothed()

    #exporting to PDF DXF and PNG
    export_to_image_pdf(interval=args.interval)
    create_map(maptype="isolinien",scale=scale, interval=args.interval)


    # Roof height
    # Enable exceptions
    ogr.UseExceptions()

    # Download the Swissbuilding3D data
    download_swissbuilding3d(buffered_gpkg)

    #extract roof height
    process_roof()

    # Create a map with height information
    create_map(maptype="heightinfo",scale=scale, interval=args.interval)

    print("All Done")


