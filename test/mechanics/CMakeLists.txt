# $Id$

#test files for mechanics module
IF(ENABLE_PYTHON_TESTS)
  ADD_TEST(p::mechanics::Brick8N                           ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Brick8N.py ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_CURRENT_SOURCE_DIR})
  ADD_TEST(p::mechanics::ConstitutiveLaws                  ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ConstitutiveLaws.py ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_CURRENT_SOURCE_DIR})
  ADD_TEST(p::mechanics::CreateNodesAndGroupOperations     ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/CreateNodesAndGroupOperations.py ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_CURRENT_SOURCE_DIR})
  ADD_TEST(p::mechanics::ConstraintsLagrange1D             ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ConstraintsLagrange1D.py ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_CURRENT_SOURCE_DIR})
  ADD_TEST(p::mechanics::GlobalDofs                        ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/GlobalDofs.py ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_CURRENT_SOURCE_DIR})
  ADD_TEST(p::mechanics::NonlocalDamagePlasticity          ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/NonlocalDamagePlasticity.py ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_CURRENT_SOURCE_DIR})
  ADD_TEST(p::mechanics::Plane2D3N                         ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Plane2D3N.py ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_CURRENT_SOURCE_DIR})
  ADD_TEST(p::mechanics::Plane2D4N                         ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Plane2D4N.py ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_CURRENT_SOURCE_DIR})
  ADD_TEST(p::mechanics::Plane2D4N_planeStress             ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Plane2D4N_planeStress.py ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_CURRENT_SOURCE_DIR})
  ADD_TEST(p::mechanics::Section                           ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Sections.py ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_CURRENT_SOURCE_DIR})
  ADD_TEST(p::mechanics::Tetrahedron10N                    ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Tetrahedron10N.py ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_CURRENT_SOURCE_DIR})
  ADD_TEST(p::mechanics::Truss1D2N                         ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Truss1D2N.py ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_CURRENT_SOURCE_DIR})
  ADD_TEST(p::mechanics::Truss1D3N                         ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Truss1D3N.py ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_CURRENT_SOURCE_DIR})
 #ADD_TEST(p::mechanics::Tetrahedron4N                     ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Tetrahedron4N.py ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_CURRENT_SOURCE_DIR})
  ADD_TEST(p::mechanics::DeleteElemNode                    ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/DeleteElemNode.py ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_CURRENT_SOURCE_DIR})

  # set PYTHONPATH environment for all examples

  SET(MYTESTS 
  p::mechanics::Brick8N
  p::mechanics::ConstitutiveLaws
  p::mechanics::CreateNodesAndGroupOperations
  p::mechanics::ConstraintsLagrange1D
  p::mechanics::GlobalDofs
  p::mechanics::NonlocalDamagePlasticity
  p::mechanics::Plane2D3N
  p::mechanics::Plane2D4N
  p::mechanics::Plane2D4N_planeStress
  p::mechanics::Section
  p::mechanics::Tetrahedron10N
  p::mechanics::Truss1D2N
  p::mechanics::Truss1D3N
  #p::mechanics::Tetrahedron4N
  p::mechanics::DeleteElemNode
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

  # GradientDamagePlasticity
  ADD_EXECUTABLE(GradientDamagePlasticity GradientDamagePlasticity.cpp)
  TARGET_LINK_LIBRARIES(GradientDamagePlasticity NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} ${ANN_LIBRARIES})
  IF (ENABLE_VISUALIZE)
    TARGET_LINK_LIBRARIES(GradientDamagePlasticity NuToVisualize)
  ENDIF (ENABLE_VISUALIZE)
  IF(MUMPS_FOUND)
    TARGET_LINK_LIBRARIES(GradientDamagePlasticity ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
  ENDIF(MUMPS_FOUND)
  IF(MINGW)
    SET_TARGET_PROPERTIES(GradientDamagePlasticity PROPERTIES LINK_FLAGS -Wl,--enable-auto-import)
  ENDIF(MINGW)
  ADD_TEST(c::mechanics::GradientDamagePlasticity ${CMAKE_CURRENT_BUILD_DIR}/GradientDamagePlasticity)

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

  # test triangular 2D 3 node element
  ADD_EXECUTABLE(Plane2D3N Plane2D3N.cpp)
  TARGET_LINK_LIBRARIES(Plane2D3N NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} )
  IF (ENABLE_VISUALIZE)
    TARGET_LINK_LIBRARIES(Plane2D3N NuToVisualize)
  ENDIF (ENABLE_VISUALIZE)
  IF(MUMPS_FOUND)
    TARGET_LINK_LIBRARIES(Plane2D3N ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
  ENDIF(MUMPS_FOUND)
  ADD_TEST(c::mechanics::Plane2D3N ${CMAKE_CURRENT_BUILD_DIR}/Plane2D3N)

  # test surface load integration in 3D
  ADD_EXECUTABLE(SurfaceLoad2D SurfaceLoad2D.cpp)
  TARGET_LINK_LIBRARIES(SurfaceLoad2D NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} )
  IF (ENABLE_VISUALIZE)
    TARGET_LINK_LIBRARIES(SurfaceLoad2D NuToVisualize)
  ENDIF (ENABLE_VISUALIZE)
  IF(MUMPS_FOUND)
    TARGET_LINK_LIBRARIES(SurfaceLoad2D ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
  ENDIF(MUMPS_FOUND)
  ADD_TEST(c::mechanics::SurfaceLoad2D ${CMAKE_CURRENT_BUILD_DIR}/SurfaceLoad2D)

  # test triangular 10 node element
  ADD_EXECUTABLE(Plane2D10N Plane2D10N.cpp)
  TARGET_LINK_LIBRARIES(Plane2D10N NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} )
  IF (ENABLE_VISUALIZE)
    TARGET_LINK_LIBRARIES(Plane2D10N NuToVisualize)
  ENDIF (ENABLE_VISUALIZE)
  IF(MUMPS_FOUND)
    TARGET_LINK_LIBRARIES(Plane2D10N ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
  ENDIF(MUMPS_FOUND)
  IF(MINGW)
    SET_TARGET_PROPERTIES(Plane2D10N PROPERTIES LINK_FLAGS -Wl,--enable-auto-import)
  ENDIF(MINGW)#
  ADD_TEST(c::mechanics::Plane2D10N ${CMAKE_CURRENT_BUILD_DIR}/Plane2D10N)
  
  # test triangular 15 node element
  ADD_EXECUTABLE(Plane2D15N Plane2D15N.cpp)
  TARGET_LINK_LIBRARIES(Plane2D15N NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} )
  IF (ENABLE_VISUALIZE)
    TARGET_LINK_LIBRARIES(Plane2D15N NuToVisualize)
  ENDIF (ENABLE_VISUALIZE)
  IF(MUMPS_FOUND)
    TARGET_LINK_LIBRARIES(Plane2D15N ${MUMPS_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
  ENDIF(MUMPS_FOUND)
  ADD_TEST(c::mechanics::Plane2D15N ${CMAKE_CURRENT_BUILD_DIR}/Plane2D15N)

  # test the CheckStiffness method globally and locally
  ADD_EXECUTABLE(CoefficientChecks CoefficientChecks.cpp)
  TARGET_LINK_LIBRARIES(CoefficientChecks NuToMechanics NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} )
  ADD_TEST(c::mechanics::CoefficientChecks ${CMAKE_CURRENT_BUILD_DIR}/CoefficientChecks)

  # ooIntegration
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
  ENDIF(MINGW)#
  ADD_TEST(c::mechanics::LobattoIntegration ${CMAKE_CURRENT_BUILD_DIR}/LobattoIntegration)

  # MoistureTransport
  IF(ENABLE_METAMODEL)
    ADD_EXECUTABLE(MoistureTransport1D MoistureTransport.cpp)
    TARGET_LINK_LIBRARIES(MoistureTransport1D NuToMetamodel NuToMechanics NuToMath NuToBase)
    ADD_TEST(c::mechanics::MoistureTransport1D ${CMAKE_CURRENT_BUILD_DIR}/MoistureTransport1D)
  ENDIF(ENABLE_METAMODEL)

  
ENDIF(ENABLE_C++_TESTS)
  
IF(ENABLE_C++_TESTS_WIP)  
  
ENDIF(ENABLE_C++_TESTS_WIP)
