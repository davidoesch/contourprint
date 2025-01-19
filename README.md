# ContourPrint

This repository contains a Python script to process geospatial data, generate contour lines, smooth geometries, and extract roof heights. It also creates maps with height information and exports the results to various formats including PDF, DXF, and PNG.

## Workflow

1. Buffer the area of interest (AOI) by 10 meters.
2. Clip the data using the buffered area.
3. Generate contour lines.
4. Smooth the contours with 10 iterations and an offset of 0.125.
5. Apply styles to the smoothed contours.

For more details, refer to:
- [Creating Contour Lines and Labels with QGIS](https://opensourceoptions.com/how-to-create-contour-lines-and-labels-with-qgis/)
- [Copying Styles in QGIS](https://gis.stackexchange.com/a/301965/190185)

## Example Call

Run the script from the command line as follows:
```bash
.venv/bin/python create_rectangle_geopackage.py --bbox "2'604’333/1'197'467 2'604’543/1'197'467 2'604’543/1'197'317 2'604’333/1'197'317" --interval 0.20

Arguments
--bbox (str): Bounding box coordinates in the format X/Y, separated by spaces. Default: "2'604’333/1'197'467 2'604’543/1'197'467 2'604’543/1'197'317 2'604’333/1'197'317"
--interval (float): Interval in meters for contour lines (äquidisztanz). Default is 0.2m.
--resampling (float): Resampling in meters. Default is 10.0.

## Requirements

Python 3.x
QGIS must be installed and accessible from the command line.
Required Python packages (install via requirements.txt):
    shapely
    pyproj
    rasterio
    matplotlib_scalebar
    osgeo
    owslib
    PIL (Pillow)
    geopandas
    pandas
    numpy

