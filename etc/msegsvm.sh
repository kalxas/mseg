./msegcli_svm toposys2.ers output WeightsList.txt LevelParam.txt
gdal_polygonize.py raster.ers -f "ESRI Shapefile" polygons.shp
mv stats.test toposys2.test
mv training.train toposys2.train
./easy.py toposys2.train toposys2.test
rm toposys2.test.scale
rm toposys2.train.range
rm toposys2.train.scale
rm toposys2.train.scale.out
./obia-class-export-svm raster TTAMask toposys2.test toposys2.test.predict classification classraster
gdal_polygonize.py classraster.ers -f "ESRI Shapefile" classes.shp
./obia-utils-confmatrix classraster truth 4 474 473 > confmatrix.txt