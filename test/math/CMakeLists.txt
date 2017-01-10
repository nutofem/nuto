if(ENABLE_PYTHON)
    foreach(test
            SparseMatrixCSR
            Numpy
            )

        add_test(p::math::${test} ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${test}.py ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR} ${CMAKE_CURRENT_SOURCE_DIR})
        set_property(TEST p::math::${test} PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python")
    endforeach()
endif()

if(ENABLE_C++_TESTS)
    add_executable(SparseMatrixTests SparseMatrixTests.cpp)
    target_link_libraries(SparseMatrixTests NuToMath NuToBase ${Boost_LIBRARIES})
    add_test(c::math::SparseMatrixTests ${CMAKE_CURRENT_BUILD_DIR}/SparseMatrixTests)

    add_executable(Interpolation Interpolation.cpp)
    target_link_libraries(Interpolation NuToMath NuToBase)
    add_test(c::unit::Interpolation ${CMAKE_CURRENT_BUILD_DIR}/Interpolation)

    add_executable(SerializeData SerializeData.cpp)
    target_link_libraries(SerializeData NuToMath NuToBase)
    add_test(c::unit::SerializeData ${CMAKE_CURRENT_BUILD_DIR}/SerializeData)

    if(ENABLE_ARPACK)
        add_executable(Arpack Arpack.cpp)
        target_link_libraries(Arpack NuToMath NuToBase)
        add_test(c::math::Arpack ${CMAKE_CURRENT_BUILD_DIR}/Arpack)
    endif()

    add_executable(EigenCompanion EigenCompanion.cpp ${CMAKE_SOURCE_DIR}/src/math/EigenCompanion.cpp)
    target_link_libraries(EigenCompanion ${Boost_LIBRARIES})
    add_test(c::unit::EigenCompanion ${CMAKE_CURRENT_BUILD_DIR}/EigenCompanion)

endif()
