# $Id$

#test files for mechanics module
IF(ENABLE_PYTHON_TESTS)
  ADD_TEST(p::mechanics::Brick8N                           ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Brick8N.py ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_CURRENT_SOURCE_DIR})
  ADD_TEST(p::mechanics::ConstitutiveLaws                  ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ConstitutiveLaws.py ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_CURRENT_SOURCE_DIR})
  ADD_TEST(p::mechanics::CreateNodesAndGroupOperations     ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/CreateNodesAndGroupOperations.py ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_CURRENT_SOURCE_DIR})
  ADD_TEST(p::mechanics::ConstraintsLagrange1D             ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ConstraintsLagrange1D.py ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_CURRENT_SOURCE_DIR})
  ADD_TEST(p::mechanics::GlobalDofs                        ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/GlobalDofs.py ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_CURRENT_SOURCE_DIR})
  ADD_TEST(p::mechanics::NonlocalDamagePlasticity          ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/NonlocalDamagePlasticity.py ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_CURRENT_SOURCE_DIR})
  ADD_TEST(p::mechanics::Section                           ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Sections.py ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_CURRENT_SOURCE_DIR})
  ADD_TEST(p::mechanics::Tetrahedron10N                    ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Tetrahedron10N.py ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_CURRENT_SOURCE_DIR})
  ADD_TEST(p::mechanics::DeleteElemNode                    ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/DeleteElemNode.py ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_CURRENT_SOURCE_DIR})
  ADD_TEST(p::mechanics::ElementUniaxialTest               ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ElementUniaxialTest.py ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_CURRENT_BINARY_DIR})
  ADD_TEST(p::mechanics::Quad2DPatchTest                   ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Quad2DPatchTest.py ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_CURRENT_BINARY_DIR})

  # set PYTHONPATH environment for all examples

  SET(MYTESTS 
  p::mechanics::Brick8N
  p::mechanics::ConstitutiveLaws
  p::mechanics::CreateNodesAndGroupOperations
  p::mechanics::ConstraintsLagrange1D
  p::mechanics::GlobalDofs
  p::mechanics::NonlocalDamagePlasticity
  p::mechanics::Section
  p::mechanics::Tetrahedron10N
  p::mechanics::DeleteElemNode
  p::mechanics::ElementUniaxialTest
  p::mechanics::Quad2DPatchTest
  )
  SET_PROPERTY(TEST ${MYTESTS}  PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/src")
ENDIF(ENABLE_PYTHON_TESTS)


IF(ENABLE_C++_TESTS)
  # Test ConstraintsLagrange
  ADD_EXECUTABLE(ConstraintsLagrange2D ConstraintsLagrange2D.cpp)
  TARGET_LINK_LIBRARIES(ConstraintsLagrange2D NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${BLAS_LIBRARIES} ${ANN_LIBRARIES})
  IF (ENABLE_VISUALIZE)
    TARGET_LINK_LIBRARIES(ConstraintsLagrange2D NuToVisualize)
  ENDIF (ENABLE_VISUALIZE)
  IF(MUMPS_FOUND)
    TARGET_LINK_LIBRARIES(ConstraintsLagrange2D ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
  ENDIF(MUMPS_FOUND)
  IF(MINGW)
    SET_TARGET_PROPERTIES(ConstraintsLagrange2D PROPERTIES LINK_FLAGS -Wl,--enable-auto-import)
  ENDIF(MINGW)
  ADD_TEST(c::mechanics::ConstraintsLagrange2D ${CMAKE_CURRENT_BUILD_DIR}/ConstraintsLagrange2D)

# GradientDamage
  ADD_EXECUTABLE(GradientDamage GradientDamage.cpp)
  TARGET_LINK_LIBRARIES(GradientDamage NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
  IF (ENABLE_VISUALIZE)
    TARGET_LINK_LIBRARIES(GradientDamage NuToVisualize)
  ENDIF (ENABLE_VISUALIZE)
  IF(MUMPS_FOUND)
    TARGET_LINK_LIBRARIES(GradientDamage ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
  ENDIF(MUMPS_FOUND)
  IF(MINGW)
    SET_TARGET_PROPERTIES(GradientDamage PROPERTIES LINK_FLAGS -Wl,--enable-auto-import)
  ENDIF(MINGW)
  ADD_TEST(c::mechanics::GradientDamage ${CMAKE_CURRENT_BUILD_DIR}/GradientDamage)

  # test Newmark time integration
  ADD_EXECUTABLE(NewmarkPlane2D4N NewmarkPlane2D4N.cpp)
  TARGET_LINK_LIBRARIES(NewmarkPlane2D4N NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
  IF (ENABLE_VISUALIZE)
    TARGET_LINK_LIBRARIES(NewmarkPlane2D4N NuToVisualize)
  ENDIF (ENABLE_VISUALIZE)
  IF(MUMPS_FOUND)
    TARGET_LINK_LIBRARIES(NewmarkPlane2D4N ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
  ENDIF(MUMPS_FOUND)
  IF(MINGW)
    SET_TARGET_PROPERTIES(NewmarkPlane2D4N PROPERTIES LINK_FLAGS -Wl,--enable-auto-import)
  ENDIF(MINGW)
  ADD_TEST(c::mechanics::NewmarkPlane2D4N ${CMAKE_CURRENT_BUILD_DIR}/NewmarkPlane2D4N)

  # test Newmark time integration with a dogbone specimen
  ADD_EXECUTABLE(DogBoneDirectNewmark DogBoneDirectNewmark.cpp)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/DogBoneDirectNewmarkGeometry.geo ${CMAKE_CURRENT_BINARY_DIR}/DogBoneDirectNewmarkGeometry.geo  COPYONLY)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/DogBoneDirectNewmarkGeometry.msh ${CMAKE_CURRENT_BINARY_DIR}/DogBoneDirectNewmarkGeometry.msh  COPYONLY)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/results/Linuxx86_64/DogboneDirectNewmarkresultAllLoadStepsRef.dat ${CMAKE_CURRENT_BINARY_DIR}/DogboneDirectNewmarkresultAllLoadStepsRef.dat  COPYONLY)
  
  TARGET_LINK_LIBRARIES(DogBoneDirectNewmark NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
  IF (ENABLE_VISUALIZE)
    TARGET_LINK_LIBRARIES(DogBoneDirectNewmark NuToVisualize)
  ENDIF (ENABLE_VISUALIZE)
  IF(MUMPS_FOUND)
    TARGET_LINK_LIBRARIES(DogBoneDirectNewmark ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
  ENDIF(MUMPS_FOUND)
  IF(MINGW)
    SET_TARGET_PROPERTIES(DogBoneDirectNewmark PROPERTIES LINK_FLAGS -Wl,--enable-auto-import)
  ENDIF(MINGW)
  ADD_TEST(c::mechanics::DogBoneDirectNewmark ${CMAKE_CURRENT_BUILD_DIR}/DogBoneDirectNewmark)

  # test surface load integration in 2D
  ADD_EXECUTABLE(SurfaceLoad SurfaceLoad.cpp)
  TARGET_LINK_LIBRARIES(SurfaceLoad NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} )
  IF (ENABLE_VISUALIZE)
    TARGET_LINK_LIBRARIES(SurfaceLoad NuToVisualize)
  ENDIF (ENABLE_VISUALIZE)
  IF(MUMPS_FOUND)
    TARGET_LINK_LIBRARIES(SurfaceLoad ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
  ENDIF(MUMPS_FOUND)
  ADD_TEST(c::mechanics::SurfaceLoad ${CMAKE_CURRENT_BUILD_DIR}/SurfaceLoad)

  # test the CheckStiffness method globally and locally
  ADD_EXECUTABLE(CoefficientChecks CoefficientChecks.cpp)
  TARGET_LINK_LIBRARIES(CoefficientChecks NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} )
  ADD_TEST(c::mechanics::CoefficientChecks ${CMAKE_CURRENT_BUILD_DIR}/CoefficientChecks)

  ADD_EXECUTABLE(ElementUniaxialTruss1D ElementUniaxialTruss1D.cpp)
  TARGET_LINK_LIBRARIES(ElementUniaxialTruss1D NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} )
  IF (ENABLE_VISUALIZE)
    TARGET_LINK_LIBRARIES(ElementUniaxialTruss1D NuToVisualize)
  ENDIF (ENABLE_VISUALIZE)
  IF(MUMPS_FOUND)
    TARGET_LINK_LIBRARIES(ElementUniaxialTruss1D ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
  ENDIF(MUMPS_FOUND)
  IF(MINGW)
    SET_TARGET_PROPERTIES(ElementUniaxialTruss1D PROPERTIES LINK_FLAGS -Wl,--enable-auto-import)
  ENDIF(MINGW)
  ADD_TEST(c::mechanics::ElementUniaxialTruss1D ${CMAKE_CURRENT_BUILD_DIR}/ElementUniaxialTruss1D)

  ADD_EXECUTABLE(ElementUniaxialTriangle2D ElementUniaxialTriangle2D.cpp)
  TARGET_LINK_LIBRARIES(ElementUniaxialTriangle2D NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} )
  IF (ENABLE_VISUALIZE)
    TARGET_LINK_LIBRARIES(ElementUniaxialTriangle2D NuToVisualize)
  ENDIF (ENABLE_VISUALIZE)
  IF(MUMPS_FOUND)
    TARGET_LINK_LIBRARIES(ElementUniaxialTriangle2D ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
  ENDIF(MUMPS_FOUND)
  IF(MINGW)
    SET_TARGET_PROPERTIES(ElementUniaxialTriangle2D PROPERTIES LINK_FLAGS -Wl,--enable-auto-import)
  ENDIF(MINGW)
  ADD_TEST(c::mechanics::ElementUniaxialTriangle2D ${CMAKE_CURRENT_BUILD_DIR}/ElementUniaxialTriangle2D)

  ADD_EXECUTABLE(ElementUniaxialQuad2D ElementUniaxialQuad2D.cpp)
  TARGET_LINK_LIBRARIES(ElementUniaxialQuad2D NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} )
  IF (ENABLE_VISUALIZE)
    TARGET_LINK_LIBRARIES(ElementUniaxialQuad2D NuToVisualize)
  ENDIF (ENABLE_VISUALIZE)
  IF(MUMPS_FOUND)
    TARGET_LINK_LIBRARIES(ElementUniaxialQuad2D ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
  ENDIF(MUMPS_FOUND)
  IF(MINGW)
    SET_TARGET_PROPERTIES(ElementUniaxialQuad2D PROPERTIES LINK_FLAGS -Wl,--enable-auto-import)
  ENDIF(MINGW)
  ADD_TEST(c::mechanics::ElementUniaxialQuad2D ${CMAKE_CURRENT_BUILD_DIR}/ElementUniaxialQuad2D)

  ADD_EXECUTABLE(ElementUniaxialTetrahedron3D ElementUniaxialTetrahedron3D.cpp)
  TARGET_LINK_LIBRARIES(ElementUniaxialTetrahedron3D NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} )
  IF (ENABLE_VISUALIZE)
    TARGET_LINK_LIBRARIES(ElementUniaxialTetrahedron3D NuToVisualize)
  ENDIF (ENABLE_VISUALIZE)
  IF(MUMPS_FOUND)
    TARGET_LINK_LIBRARIES(ElementUniaxialTetrahedron3D ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
  ENDIF(MUMPS_FOUND)
  IF(MINGW)
    SET_TARGET_PROPERTIES(ElementUniaxialTetrahedron3D PROPERTIES LINK_FLAGS -Wl,--enable-auto-import)
  ENDIF(MINGW)
  ADD_TEST(c::mechanics::ElementUniaxialTetrahedron3D ${CMAKE_CURRENT_BUILD_DIR}/ElementUniaxialTetrahedron3D)

  ADD_EXECUTABLE(ElementUniaxialBrick3D ElementUniaxialBrick3D.cpp)
  TARGET_LINK_LIBRARIES(ElementUniaxialBrick3D NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} )
  IF (ENABLE_VISUALIZE)
    TARGET_LINK_LIBRARIES(ElementUniaxialBrick3D NuToVisualize)
  ENDIF (ENABLE_VISUALIZE)
  IF(MUMPS_FOUND)
    TARGET_LINK_LIBRARIES(ElementUniaxialBrick3D ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
  ENDIF(MUMPS_FOUND)
  IF(MINGW)
    SET_TARGET_PROPERTIES(ElementUniaxialBrick3D PROPERTIES LINK_FLAGS -Wl,--enable-auto-import)
  ENDIF(MINGW)
  ADD_TEST(c::mechanics::ElementUniaxialBrick3D ${CMAKE_CURRENT_BUILD_DIR}/ElementUniaxialBrick3D)

  ADD_EXECUTABLE(LobattoIntegration LobattoIntegration.cpp)
  TARGET_LINK_LIBRARIES(LobattoIntegration NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} )
  IF (ENABLE_VISUALIZE)
    TARGET_LINK_LIBRARIES(LobattoIntegration NuToVisualize)
  ENDIF (ENABLE_VISUALIZE)
  IF(MUMPS_FOUND)
    TARGET_LINK_LIBRARIES(LobattoIntegration ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
  ENDIF(MUMPS_FOUND)
  IF(MINGW)
    SET_TARGET_PROPERTIES(LobattoIntegration PROPERTIES LINK_FLAGS -Wl,--enable-auto-import)
  ENDIF(MINGW)
  ADD_TEST(c::mechanics::LobattoIntegration ${CMAKE_CURRENT_BUILD_DIR}/LobattoIntegration)

# NonlocalDamagePlasticityModel
ADD_EXECUTABLE(NonlocalDamagePlasticityModel NonlocalDamagePlasticityModel.cpp)
TARGET_LINK_LIBRARIES(NonlocalDamagePlasticityModel NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
IF (ENABLE_VISUALIZE)
  TARGET_LINK_LIBRARIES(NonlocalDamagePlasticityModel NuToVisualize)
ENDIF (ENABLE_VISUALIZE)
IF(MUMPS_FOUND)
  TARGET_LINK_LIBRARIES(NonlocalDamagePlasticityModel ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
ENDIF(MUMPS_FOUND)
IF(MINGW)
  SET_TARGET_PROPERTIES(NonlocalDamagePlasticityModel PROPERTIES LINK_FLAGS -Wl,--enable-auto-import)
ENDIF(MINGW)
ADD_TEST(c::mechanics::NonlocalDamagePlasticityModel ${CMAKE_CURRENT_BUILD_DIR}/NonlocalDamagePlasticityModel)

  # MoistureTransport
  IF(ENABLE_METAMODEL)
    ADD_EXECUTABLE(MoistureTransport1D MoistureTransport1D.cpp)
    TARGET_LINK_LIBRARIES(MoistureTransport1D NuToMetamodel NuToMechanics NuToMath NuToBase)
    ADD_TEST(c::mechanics::MoistureTransport1D ${CMAKE_CURRENT_BUILD_DIR}/MoistureTransport1D)
  IF(MINGW)
    SET_TARGET_PROPERTIES(MoistureTransport1D PROPERTIES LINK_FLAGS -Wl,--enable-auto-import)
  ENDIF(MINGW)
  ENDIF(ENABLE_METAMODEL)

  #  new interpolation type tests
  CONFIGURE_FILE(2D.msh ${CMAKE_CURRENT_BINARY_DIR})
  ADD_EXECUTABLE(InterpolationTypes InterpolationTypes.cpp)
  TARGET_LINK_LIBRARIES(InterpolationTypes NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} )
  IF (ENABLE_VISUALIZE)
    TARGET_LINK_LIBRARIES(InterpolationTypes NuToVisualize)
  ENDIF (ENABLE_VISUALIZE)
  IF(MUMPS_FOUND)
    TARGET_LINK_LIBRARIES(InterpolationTypes ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
  ENDIF(MUMPS_FOUND)
  ADD_TEST(c::mechanics::InterpolationTypes ${CMAKE_CURRENT_BUILD_DIR}/InterpolationTypes)


    # structure evaluate
  ADD_EXECUTABLE(StructureEvaluate StructureEvaluate.cpp)
    TARGET_LINK_LIBRARIES(StructureEvaluate NuToMechanics NuToMath NuToBase)
    ADD_TEST(c::mechanics::StructureEvaluate ${CMAKE_CURRENT_BUILD_DIR}/StructureEvaluate)
  IF(MINGW)
    SET_TARGET_PROPERTIES(StructureEvaluate PROPERTIES LINK_FLAGS -Wl,--enable-auto-import)
  ENDIF(MINGW)
  
ENDIF(ENABLE_C++_TESTS)
  
IF(ENABLE_C++_TESTS_WIP)  
 
ENDIF(ENABLE_C++_TESTS_WIP)
