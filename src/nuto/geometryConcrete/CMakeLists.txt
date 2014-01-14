# $Id: CMakeLists.txt 657 2014-01-13 09:34:14Z unger3 $
SET(NuToGeometryConcreteCollision_SOURCES
collision/HelloWorld.cpp
)

SET(NuToGeometryConcrete_SOURCES
${NuToGeometryConcreteCollision_SOURCES}
)

# build static library
ADD_LIBRARY(NuToGeometryConcreteCollision ${NuToGeometryConcreteCollision_SOURCES})
TARGET_LINK_LIBRARIES(NuToGeometryConcreteCollision ${Boost_LIBRARIES})
IF(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
  SET_TARGET_PROPERTIES(NuToGeometryConcreteCollision PROPERTIES COMPILE_FLAGS "-fPIC")
ENDIF(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )

# install library
INSTALL(TARGETS NuToGeometryConcreteCollision
        ARCHIVE DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH}
	LIBRARY DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH})

# build python modules
#CONFIGURE_FILE(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py)
#NUTO_INSTALL_PYTHON_FILE(__init__.py ${NUTO_PYTHON_MODULES_INSTALL_PATH}/nuto/base)
#FOREACH(NuToBase_SWIG_MODULE ${NuToBase_SWIG_MODULES})
#  NUTO_SWIG_MODULE(${NuToBase_SWIG_MODULE} ${NuToBase_SWIG_MODULE}.i nuto/base NuToBase)
#ENDFOREACH(NuToBase_SWIG_MODULE)

