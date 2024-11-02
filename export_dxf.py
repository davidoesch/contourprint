import os
from qgis.core import QgsApplication, QgsProject, QgsVectorLayer, QgsVectorFileWriter, QgsCoordinateReferenceSystem, QgsRectangle

# Suppress GUI-related warnings
os.environ['SESSION_MANAGER'] = ''
os.environ['QT_QPA_PLATFORM'] = 'offscreen'

# Initialize QGIS application
app = QgsApplication([], False)
app.initQgis()

# Create a new project
project = QgsProject.instance()

# Set the input and output paths
input_gpkg = "/media/menas/data/projects/contourprint/swissalti3d_data/contour_smoothed.gpkg"
output_dxf = "/media/menas/data/projects/contourprint/swissalti3d_data/test_output.dxf"

# Load the layer from the GeoPackage
layer = QgsVectorLayer(f"{input_gpkg}|layername=contour_smoothed", "Contours", "ogr")

if not layer.isValid():
    print("Layer failed to load!")
else:
    # Add layer to the project
    project.addMapLayer(layer)
    print("Layer loaded successfully")

    # Set up the DXF export options
    options = QgsVectorFileWriter.SaveVectorOptions()
    options.driverName = "DXF"
    options.fileEncoding = "ISO-8859-1"
    options.ct = QgsCoordinateReferenceSystem("EPSG:2056")
    options.layerOptions = ['FORCE2D=NO', 'MTEXT=YES']
    options.symbologyExport = QgsVectorFileWriter.NoSymbology
    options.symbologyScale = 500.0
    options.filterExtent = None  # Set to QgsRectangle() if you want to specify an extent

    # Export to DXF
    error = QgsVectorFileWriter.writeAsVectorFormat(layer, output_dxf, options)

    if error[0] == QgsVectorFileWriter.NoError:
        print(f"Layer exported successfully to {output_dxf}")
    else:
        print(f"Error exporting layer: {error}")

# Clean up
app.exitQgis()