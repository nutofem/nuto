#test files for math module
IF(ENABLE_PYTHON_TESTS)
IF(ENABLE_METAMODEL)
  SET( CMAKE_VERBOSE_MAKEFILE on )

  ADD_TEST(p::metamodel::NeuralNetwork ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/NeuralNetwork.py ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_CURRENT_SOURCE_DIR})
  SET_PROPERTY(TEST p::metamodel::NeuralNetwork PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/src")
ENDIF(ENABLE_METAMODEL)
ENDIF(ENABLE_PYTHON_TESTS)

IF(ENABLE_C++_TESTS)
IF(ENABLE_METAMODEL)
  # test the CheckStiffness method globally and locally
  ADD_EXECUTABLE(PolynomialLeastSquaresFitting PolynomialLeastSquaresFitting.cpp)
  TARGET_LINK_LIBRARIES(PolynomialLeastSquaresFitting NuToMetamodel NuToMath NuToBase ${Boost_LIBRARIES} ${LAPACK_LIBRARIES} )
  ADD_TEST(c::metamodel::PolynomialLeastSquaresFitting ${CMAKE_CURRENT_BUILD_DIR}/PolynomialLeastSquaresFitting)
ENDIF(ENABLE_METAMODEL)
ENDIF(ENABLE_C++_TESTS)
