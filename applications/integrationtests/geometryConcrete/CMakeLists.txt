if(ENABLE_GEOMETRYCONCRETE)
    unset(GMSH_PATH CACHE)
    find_program(GMSH_PATH gmsh)
    if(ENABLE_PYTHON)
        foreach(inputFile
            input2D_test.dat
            input3D_test.dat
            )
            create_symlink("${inputFile}")
        endforeach()

        foreach(test
            CollisionDetection
            GradingCurveFileIO
            )

            add_test(integration::geometryConcrete::${test}.py ${CMAKE_CURRENT_SOURCE_DIR}/${test}.py)
            set_property(TEST integration::geometryConcrete::${test}.py PROPERTY ENVIRONMENT
                "PYTHONPATH=${CMAKE_BINARY_DIR}/python")
        endforeach()

        if(GMSH_PATH)
            if(ENABLE_MECHANICS)
                add_test(integration::geometryConcrete::GmshExport.py ${CMAKE_CURRENT_SOURCE_DIR}/GmshExport.py
                    ${GMSH_PATH})
                set_property(TEST integration::geometryConcrete::GmshExport.py PROPERTY ENVIRONMENT
                    "PYTHONPATH=${CMAKE_BINARY_DIR}/python")
            endif()

        else()
            message(STATUS "Gmsh not found. Skipping test integration::geometryConcrete::GmshExport.py")
        endif()
    endif()

    if(GMSH_PATH)
        add_executable(GeometryConcreteAPI GeometryConcreteAPI.cpp)
        target_link_libraries(GeometryConcreteAPI NuToMath NuToBase NuToGeometryConcrete ${Boost_LIBRARIES})
        add_test(integration::geometryConcrete::GeometryConcreteAPI.cpp
            ${CMAKE_CURRENT_BUILD_DIR}/GeometryConcreteAPI ${GMSH_PATH})
    else()
        message(STATUS "Gmsh not found. Skipping test integration::geometryConcrete::GeometryConcreteAPI.cpp")
    endif()
endif()
