Welcome to Mseg v0.9.5

Mseg is a generic region-based multi-scale image segmentation algorithm designed with 
some optimizations for remote sensing imagery. 

The algorithm can be used as a low level processing part of an free object-oriented image analysis system. 

Mseg is Free Software licenced under the GPLv2.

Installation:

1. Dependencies: CMake, FreeImage, TinyXML

2. Run cmake configuration:

cmake -DCMAKE_BUILD_TYPE:STRING="Release" \		#or "Debug"
      -DCMAKE_INSTALL_PREFIX:PATH=/usr \
      -DLIB_SUFFIX="64" 				#use in case of 64bit OS

3. Compile the sources:

make

4. Install:

make install

For problems / suggestions / bugs etc please contact Angelos Tzotsos (tzotsos@gmail.com)