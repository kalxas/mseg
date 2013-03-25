#
# Try to find the FreeImage library and include path.
# Once done this will define
#
# FREEIMAGE_FOUND
# FREEIMAGE_INCLUDE_PATH
# FREEIMAGE_LIBRARY
# 

FIND_PATH( TINYXML_INCLUDE_PATH tinyxml.h
	/usr/include
	/usr/local/include
	/sw/include
	/opt/local/include
	DOC "The directory where TinyXML.h resides")
FIND_LIBRARY( TINYXML_LIBRARY
	NAMES TinyXML tinyxml
	PATHS
	/usr/lib64
	/usr/lib
	/usr/local/lib64
	/usr/local/lib
	/sw/lib
	/opt/local/lib
	DOC "The TinyXML library")

SET(TINYXML_LIBRARIES ${TINYXML_LIBRARY})

MARK_AS_ADVANCED(
	TINYXML_LIBRARY
	TINYXML_LIBRARIES
	TINYXML_INCLUDE_PATH)