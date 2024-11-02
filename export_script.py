import os
os.environ['QT_QPA_PLATFORM'] = 'offscreen'
print("QGIS initialized")
from qgis.core import QgsApplication, QgsProject, QgsLayoutExporter, QgsVectorLayer, QgsLineSymbol, QgsPalLayerSettings, QgsTextFormat, QgsVectorLayerSimpleLabeling, QgsProperty,QgsExpression,QgsExpressionContext, QgsExpressionContextUtils
from qgis.PyQt.QtCore import QSize
from qgis.PyQt.QtGui import QFont, QColor

app = QgsApplication([], False)
app.initQgis()


project = QgsProject.instance()
project.read('/media/menas/data/projects/contourprint/template2.qgz')
print("Project loaded")


# Load and add contour layer
contour_layer = QgsVectorLayer("/media/menas/data/projects/contourprint/swissalti3d_data/contour_smoothed.gpkg", "contour_smoothed", "ogr")
if not contour_layer.isValid():
    print("Contour layer failed to load!")
else:
    project.addMapLayer(contour_layer)
    print("Contour layer added")

    # Style the contour layer
    symbol = QgsLineSymbol.createSimple({'line_color': 'black'})

    # Set line width based on elevation
    width_expression = 'if("ELEV" % 2 = 0, 1.8, if("ELEV" % 1 = 0, 0.9, 1.36))'
    symbol.setDataDefinedWidth(QgsProperty.fromExpression(width_expression))

    contour_layer.renderer().setSymbol(symbol)

    # Label the contours
    layer_settings = QgsPalLayerSettings()
    text_format = QgsTextFormat()
    text_format.setFont(QFont("Arial", 28))
    text_format.setSize(28)
    text_format.setColor(QColor('black'))
    layer_settings.setFormat(text_format)

    # Set label content based on elevation
    layer_settings.fieldName = 'if("ELEV" % 1 = 0, "ELEV" || \' m\', \'\')'
    layer_settings.isExpression = True

    layer_settings.placement = QgsPalLayerSettings.Line
    layer_settings.enabled = True

    layer_labeling = QgsVectorLayerSimpleLabeling(layer_settings)
    contour_layer.setLabelsEnabled(True)
    contour_layer.setLabeling(layer_labeling)

    contour_layer.triggerRepaint()
    print("Contour layer styled and labeled")

    # Verify labeling is enabled and expression is set
    print(f"Labels enabled: {contour_layer.labelsEnabled()}")
    print(f"Labeling expression: {contour_layer.labeling().settings().fieldName}")

# Adjust map extent to zoom to contour layer
layout = project.layoutManager().layoutByName('a0v2')
map_item = layout.itemById('Map 1')  # Replace 'Map 1' with the actual ID of your map item

if map_item:
    extent = contour_layer.extent()
    extent.scale(1.1)  # Add 10% padding
    map_item.setExtent(extent)
    print("Map extent adjusted to contour layer")
else:
    print("Map item not found in layout")

exporter = QgsLayoutExporter(layout)
print("Layout found")

settings = QgsLayoutExporter.PdfExportSettings()
settings.dpi = 300  # Adjust DPI as needed
settings.scale = 1000 # Set your desired scale
print("Exporting to PDF...")

try:
    exporter.exportToPdf('/media/menas/data/projects/contourprint/output.pdf', settings)
    print("PDF export completed")
except Exception as e:
    print(f"Error exporting PDF: {str(e)}")

app.exitQgis()

