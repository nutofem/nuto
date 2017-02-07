set(NuToGeometryConcreteCollision_SOURCES
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

set(NuToGeometryConcreteGradingCurve_SOURCES
    takeAndPlace/ParticleCreator.cpp
    )

set(NuToGeometryConcrete_SOURCES
    GeometryConcrete.cpp
    InputReader.cpp
    Specimen.cpp
    ${NuToGeometryConcreteCollision_SOURCES}
    ${NuToGeometryConcreteGradingCurve_SOURCES}
    )

# build static library
add_library(NuToGeometryConcrete ${NuToGeometryConcrete_SOURCES})
if(ENABLE_VISUALIZE)
    target_link_libraries(NuToGeometryConcrete ${Boost_LIBRARIES} NuToBase NuToMath NuToVisualize)
else()
    target_link_libraries(NuToGeometryConcrete ${Boost_LIBRARIES} NuToBase NuToMath)
endif()

if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
    set_target_properties(NuToGeometryConcrete PROPERTIES COMPILE_FLAGS "-fPIC")
endif()

# install library
install(TARGETS NuToGeometryConcrete
    ARCHIVE DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH}
    LIBRARY DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH}
    )
