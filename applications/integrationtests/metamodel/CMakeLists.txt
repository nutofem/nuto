if(ENABLE_PYTHON)
    add_test(integration::metamodel::NeuralNetwork.py ${PYTHON_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/NeuralNetwork.py)
    set_property(TEST integration::metamodel::NeuralNetwork.py
        PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python")
endif()

# test the CheckStiffness method globally and locally
add_executable(PolynomialLeastSquaresFitting
    EXCLUDE_FROM_ALL PolynomialLeastSquaresFitting.cpp)
target_link_libraries(PolynomialLeastSquaresFitting
    Metamodel Math Base)
add_test(integration::metamodel::PolynomialLeastSquaresFitting.cpp
    ${CMAKE_CURRENT_BUILD_DIR}/PolynomialLeastSquaresFitting)
append_to_tests(PolynomialLeastSquaresFitting)
