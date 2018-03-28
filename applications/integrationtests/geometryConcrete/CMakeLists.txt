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

        add_test(integration::geometryConcrete::${test}.py
            ${CMAKE_CURRENT_SOURCE_DIR}/${test}.py)
        set_property(TEST integration::geometryConcrete::${test}.py
            PROPERTY ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}/python")
    endforeach()
endif()

if(GMSH_PATH)
    add_integrationtest(GeometryConcreteAPI)

    add_integrationtest(GmshWriter -- ${GMSH_PATH})
    target_link_libraries(GmshWriter Boost::filesystem)
else()
    message(STATUS
        "Gmsh not found. Skipping integration test GeometryConcreteAPI and GmshWriter.")
endif()
