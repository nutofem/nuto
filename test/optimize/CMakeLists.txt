if(ENABLE_PYTHON_TESTS)
    add_test(p::optimize::ConjugateGradientNonLinear ${PYTHON_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/ConjugateGradientNonLinear.py
        ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_CURRENT_SOURCE_DIR})
    set_property(TEST p::optimize::ConjugateGradientNonLinear PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/src")
endif()
