# $Id: CMakeLists.txt 657 2014-01-13 09:34:14Z unger3 $

SET(NuToGeometryConcreteCollision_SOURCES
collision/collidables/CollidableBase.cpp
collision/collidables/CollidableParticleBase.cpp
collision/collidables/CollidableParticleSphere.cpp
collision/collidables/CollidableWallBase.cpp
collision/collidables/CollidableWallPhysical.cpp
collision/collidables/CollidableWallVirtual.cpp
collision/collidables/CollidableWallCylinder.cpp
collision/Event.cpp
collision/SubBox.cpp
collision/handler/SubBoxHandler.cpp
collision/handler/CollisionHandler.cpp
collision/handler/EventListHandler.cpp
collision/handler/ParticleHandler.cpp
)


SET(NuToGeometryConcreteGradingCurve_SOURCES
takeAndPlace/ParticleCreator.cpp
)

SET(NuToGeometryConcrete_SOURCES
InputReader.cpp
WallTime.cpp
Specimen.cpp
${NuToGeometryConcreteCollision_SOURCES}
${NuToGeometryConcreteGradingCurve_SOURCES}
)

SET(NuToGeometryConcrete_SWIG_MODULES
  ModulGeometryConcrete
)

#SET(SWIG_MODULE_GeometryConcrete_EXTRA_DEPS 
#  structures/grid/MultiGridStructure.h
#)

# build static library
ADD_LIBRARY(NuToGeometryConcrete ${NuToGeometryConcrete_SOURCES})
IF(ENABLE_VISUALIZE)
TARGET_LINK_LIBRARIES(NuToGeometryConcrete ${Boost_LIBRARIES} NuToBase NuToMath NuToVisualize)
ELSE(ENABLE_VISUALIZE)
TARGET_LINK_LIBRARIES(NuToGeometryConcrete ${Boost_LIBRARIES} NuToBase NuToMath)
ENDIF(ENABLE_VISUALIZE)


IF(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
  SET_TARGET_PROPERTIES(NuToGeometryConcrete PROPERTIES COMPILE_FLAGS "-fPIC")
ENDIF(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )

# install library
INSTALL(TARGETS NuToGeometryConcrete
        ARCHIVE DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH}
	LIBRARY DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH})

# build python modules
CONFIGURE_FILE(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py)
NUTO_INSTALL_PYTHON_FILE(__init__.py ${NUTO_PYTHON_MODULES_INSTALL_PATH}/nuto/geometryConcrete)
FOREACH(NuToGeometryConcrete_SWIG_MODULE ${NuToGeometryConcrete_SWIG_MODULES})
  NUTO_SWIG_MODULE(${NuToGeometryConcrete_SWIG_MODULE} ${NuToGeometryConcrete_SWIG_MODULE}.i nuto/geometryConcrete NuToGeometryConcrete)
ENDFOREACH(NuToGeometryConcrete_SWIG_MODULE)

