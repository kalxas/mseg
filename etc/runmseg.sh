./msegcli CASI_8bit.ers output WeightsList.txt LevelParam.txt
gdal_polygonize.py raster.ers -f "ESRI Shapefile" polygons.shp
