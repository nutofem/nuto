if(ENABLE_PYTHON)
    add_test(integration::metamodel::NeuralNetwork.py ${PYTHON_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/NeuralNetwork.py)
    set_property(TEST integration::metamodel::NeuralNetwork.py
        PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python")
endif()

# test the CheckStiffness method globally and locally
add_integrationtest(PolynomialLeastSquaresFitting)
