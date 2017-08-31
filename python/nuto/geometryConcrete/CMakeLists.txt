set(geometryDir ${CMAKE_SOURCE_DIR}/src/geometryConcrete)
set(SWIG_MODULE_ModuleGeometryConcrete_EXTRA_DEPS
    ${geometryDir}/collision
    ${geometryDir}/Specimen.h
    ${geometryDir}/collision/collidables/CollidableBase.h
    ${geometryDir}/collision/collidables/CollidableParticleBase.h
    ${geometryDir}/collision/collidables/CollidableParticleSphere.h
    ${geometryDir}/collision/collidables/CollidableWallPhysical.h
    ${geometryDir}/collision/Event.h
    ${geometryDir}/collision/SubBox.h
    ${geometryDir}/collision/handler/SubBoxHandler.h
    ${geometryDir}/collision/handler/EventListHandler.h
    ${geometryDir}/collision/handler/ParticleHandler.h
    ${geometryDir}/collision/handler/CollisionHandler.h
    ${geometryDir}/takeAndPlace/ParticleCreator.h
    ${geometryDir}/InputReader.h
    ${CMAKE_SOURCE_DIR}/python/nuto/math/NuToMath.i
    )

configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py)
nuto_swig_module(ModuleGeometryConcrete ModuleGeometryConcrete.i
    GeometryConcrete)
