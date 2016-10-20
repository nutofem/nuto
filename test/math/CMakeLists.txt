if(ENABLE_PYTHON_TESTS)
    add_test(p::importNuto ${PYTHON_EXECUTABLE} -c "import nuto")
    set_property(TEST p::importNuto PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/src")

    add_test(p::math::Matrix ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Matrix.py ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_CURRENT_SOURCE_DIR})
    set_property(TEST p::math::Matrix PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/src")

    if(ENABLE_NUMPY)
        add_test(p::math::NumpyFullMatrix ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/NumpyFullMatrix.py ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_CURRENT_SOURCE_DIR})
        set_property(TEST p::math::NumpyFullMatrix PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/src")
    endif()

    add_test(p::math::SparseMatrixCSR ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/SparseMatrixCSR.py ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_CURRENT_SOURCE_DIR})
    set_property(TEST p::math::SparseMatrixCSR PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/src")
endif()

if(ENABLE_C++_TESTS)

    if(ENABLE_ARPACK)
        add_executable(Arpack Arpack.cpp)
        target_link_libraries(Arpack NuToMath NuToBase ${Boost_LIBRARIES} ${OpenBLAS_LIBRARIES})
        target_link_libraries(Arpack ${ARPACK_LIBRARIES})
        add_test(c::math::Arpack ${CMAKE_CURRENT_BUILD_DIR}/Arpack)
    endif()

    add_executable(SparseMatrixTests SparseMatrixTests.cpp)
    target_link_libraries(SparseMatrixTests NuToMath NuToBase ${Boost_LIBRARIES})
    add_test(c::math::SparseMatrixTests ${CMAKE_CURRENT_BUILD_DIR}/SparseMatrixTests)
endif()

if(ENABLE_UNIT_TESTS)
    add_executable(Interpolation Interpolation.cpp)
    target_link_libraries(Interpolation NuToMath NuToBase)
    add_test(c::unit::Interpolation ${CMAKE_CURRENT_BUILD_DIR}/Interpolation)

    add_executable(SerializeData SerializeData.cpp)
    target_link_libraries(SerializeData NuToMath NuToBase)
    add_test(c::unit::SerializeData ${CMAKE_CURRENT_BUILD_DIR}/SerializeData)
endif()
