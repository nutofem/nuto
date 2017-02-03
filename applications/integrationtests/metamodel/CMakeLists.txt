if(ENABLE_PYTHON)
    add_test(integration::metamodel::NeuralNetwork ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/NeuralNetwork.py ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_CURRENT_SOURCE_DIR})
    set_property(TEST integration::metamodel::NeuralNetwork PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python")
endif()

# test the CheckStiffness method globally and locally
add_executable(PolynomialLeastSquaresFitting PolynomialLeastSquaresFitting.cpp)
target_link_libraries(PolynomialLeastSquaresFitting NuToMetamodel NuToMath NuToBase ${Boost_LIBRARIES})
add_test(integration::metamodel::PolynomialLeastSquaresFitting ${CMAKE_CURRENT_BUILD_DIR}/PolynomialLeastSquaresFitting)
