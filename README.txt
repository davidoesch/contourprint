2671330/1169630

2671370/1169646

2671381/1169618

2671341/1169602

Workflow:
According https://opensourceoptions.com/how-to-create-contour-lines-and-labels-with-qgis/
buffer aoi 10m
clip with buffer
run contour
smooth with iteration 10, offset 0.125
copy styles https://gis.stackexchange.com/a/301965/190185


run it like This from the CLI, since QGIS will be called
example call: .venv/bin/python create_rectangle_geopackage.py --bbox "2'604’333/1'197'467 2'604’543/1'197'467 2'604’543/1'197'317 2'604’333/1'197'317" --interval 0.20

