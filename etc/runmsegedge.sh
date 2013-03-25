./msegcli_edge CASI_8bit.ers CASI_Canny.ers output WeightsList.txt LevelParam.txt
gdal_polygonize.py raster.ers -f "ESRI Shapefile" polygons.shp
