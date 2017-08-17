set(GeometryConcreteCollisionSources
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

set(GeometryConcreteGradingCurveSources
    takeAndPlace/ParticleCreator.cpp
    )

set(GeometryConcreteSources
    GeometryConcrete.cpp
    InputReader.cpp
    Specimen.cpp
    ${GeometryConcreteCollisionSources}
    ${GeometryConcreteGradingCurveSources}
    )

create_nuto_module(GeometryConcrete "${GeometryConcreteSources}")
if(ENABLE_VISUALIZE)
    target_link_libraries(GeometryConcrete ${Boost_LIBRARIES} Base Math Visualize)
else()
    target_link_libraries(GeometryConcrete ${Boost_LIBRARIES} Base Math)
endif()
