import argparse
import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon,MultiPolygon
from shapely.geometry import LineString, Point
from shapely.geometry import box
from shapely.ops import polygonize
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


# example call: .venv/bin/python create_rectangle_geopackage.py --bbox "2'604’333/1'197'467 2'604’543/1'197'467 2'604’543/1'197'317 2'604’333/1'197'317" --interval 0.20

global resample
global output_folder
output_folder = os.path.join(os.getcwd(), "swissalti3d_data")
BUFFER_M = 10
SWISSALTI3D_COLLECTION = "ch.swisstopo.swissalti3d"
STAC_API_URL = "https://data.geo.admin.ch/api/stac/v0.9"


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

    # Step 3: Define a function to set line widths based on elevation
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
        default=None,
        help="Resampling in meters . Default is None."
    )

    args = parser.parse_args()

    print(" -------------------")
    print()
    print("IMPORTANT: DO NOT Run from activated .venv, not from within VSCODE")
    print()
    print(" -------------------")




    resample=args.resampling
    output_folder = os.path.join(os.getcwd(), "swissalti3d_data",str(args.resampling))

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