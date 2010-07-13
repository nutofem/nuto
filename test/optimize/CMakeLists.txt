#test files for math module
SET( CMAKE_VERBOSE_MAKEFILE on )
ADD_TEST(optimize::ConjugateGradientNonLinear ${PYTHON_EXECUTABLE}
${CMAKE_CURRENT_SOURCE_DIR}/ConjugateGradientNonLinear.py ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_CURRENT_SOURCE_DIR})

