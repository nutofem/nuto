set(NuToGeometryConcreteCollisionSources
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

set(NuToGeometryConcreteGradingCurveSources
    takeAndPlace/ParticleCreator.cpp
    )

set(NuToGeometryConcreteSources
    GeometryConcrete.cpp
    InputReader.cpp
    Specimen.cpp
    ${NuToGeometryConcreteCollisionSources}
    ${NuToGeometryConcreteGradingCurveSources}
    )

sources_to_objects("${NuToGeometryConcreteSources}" NuToGeometryConcreteObjects)

# build static library
add_library(NuToGeometryConcrete ${NuToGeometryConcreteObjects})
if(ENABLE_VISUALIZE)
    target_link_libraries(NuToGeometryConcrete ${Boost_LIBRARIES} NuToBase NuToMath NuToVisualize)
else()
    target_link_libraries(NuToGeometryConcrete ${Boost_LIBRARIES} NuToBase NuToMath)
endif()

install(TARGETS NuToGeometryConcrete
    LIBRARY DESTINATION lib
    )
