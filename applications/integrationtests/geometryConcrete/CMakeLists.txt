if(ENABLE_GEOMETRYCONCRETE)
if(ENABLE_PYTHON)
    # add all examples
    add_test(p::geometryConcrete::CollisionDetection ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/CollisionDetection.py)
    add_test(p::geometryConcrete::GradingCurve3D     ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/GradingCurve3D.py)
    add_test(p::geometryConcrete::GradingCurveFileIO ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/GradingCurveFileIO.py ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

    # set PYTHONPATH environment for all examples
    set(MYTESTS
        p::geometryConcrete::CollisionDetection
        p::geometryConcrete::GradingCurveFileIO
        p::geometryConcrete::GradingCurve3D
    )
    set_property(TEST ${MYTESTS} PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python")

    unset(GMSH_PATH CACHE)
    find_program(GMSH_PATH gmsh)
    # If gmsh is not found, the variables value will be "GMSH_PATH-NOTFOUND", which
    # is recognized by the IF()-statement
    # The link to the mechanics library is required, since the successful import
    # to a structure is tested.
    if(GMSH_PATH)
        if(ENABLE_MECHANICS)
            message(STATUS "Gmsh found: ${GMSH_PATH}")
            add_test(p::geometryConcrete::GmshExport ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/GmshExport.py ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${GMSH_PATH})
            set_property(TEST p::geometryConcrete::GmshExport PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python")
        endif()

        add_executable(GeometryConcreteAPI GeometryConcreteAPI.cpp)
        target_link_libraries(GeometryConcreteAPI NuToMath NuToBase NuToGeometryConcrete ${Boost_LIBRARIES})
        add_test(c::geometryConcrete::GeometryConcreteAPI ${CMAKE_CURRENT_BUILD_DIR}/GeometryConcreteAPI ${GMSH_PATH})

    else()
        message(STATUS "Gmsh not found. Skipping test p::geometryConcrete::GmshExport")
        message(STATUS "Gmsh not found. Skipping test c::geometryConcrete::GeometryConcreteAPI")
    endif()
endif()
endif()
