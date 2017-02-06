if(ENABLE_PYTHON)
    add_test(integration::optimize::ConjugateGradientNonLinear ${PYTHON_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/ConjugateGradientNonLinear.py
        ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_CURRENT_SOURCE_DIR})
    set_property(TEST integration::optimize::ConjugateGradientNonLinear PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python")
endif()
