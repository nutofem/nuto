set(SWIG_MODULE_ModulGeometryConcrete_EXTRA_DEPS
        ${CMAKE_SOURCE_DIR}/src/geometryConcrete/collision
        ${CMAKE_SOURCE_DIR}/src/geometryConcrete/Specimen.h
        ${CMAKE_SOURCE_DIR}/src/geometryConcrete/collision/collidables/CollidableBase.h
        ${CMAKE_SOURCE_DIR}/src/geometryConcrete/collision/collidables/CollidableParticleBase.h
        ${CMAKE_SOURCE_DIR}/src/geometryConcrete/collision/collidables/CollidableParticleSphere.h
        ${CMAKE_SOURCE_DIR}/src/geometryConcrete/collision/collidables/CollidableWallPhysical.h
        ${CMAKE_SOURCE_DIR}/src/geometryConcrete/collision/Event.h
        ${CMAKE_SOURCE_DIR}/src/geometryConcrete/collision/SubBox.h
        ${CMAKE_SOURCE_DIR}/src/geometryConcrete/collision/handler/SubBoxHandler.h
        ${CMAKE_SOURCE_DIR}/src/geometryConcrete/collision/handler/EventListHandler.h
        ${CMAKE_SOURCE_DIR}/src/geometryConcrete/collision/handler/ParticleHandler.h
        ${CMAKE_SOURCE_DIR}/src/geometryConcrete/collision/handler/CollisionHandler.h
        ${CMAKE_SOURCE_DIR}/src/geometryConcrete/takeAndPlace/ParticleCreator.h
        ${CMAKE_SOURCE_DIR}/src/geometryConcrete/InputReader.h
        )



configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py)
nuto_install_python_file(__init__.py ${NUTO_PYTHON_MODULES_INSTALL_PATH}/nuto/geometryConcrete)
nuto_swig_module(ModulGeometryConcrete ModulGeometryConcrete.i nuto/geometryConcrete NuToGeometryConcrete)


