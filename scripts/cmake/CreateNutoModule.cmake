function(create_nuto_module ModuleName ModuleSources)
    add_library(${ModuleName} ${ModuleSources} ${ARGN})
    set_target_properties(${ModuleName} PROPERTIES
        OUTPUT_NAME NuTo${ModuleName}
        )
    target_include_directories(${ModuleName} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
        $<INSTALL_INTERFACE:include/nuto>
        )

    install(TARGETS ${ModuleName} EXPORT NuToTargets
        LIBRARY DESTINATION lib
        INCLUDES DESTINATION include/nuto
        )
endfunction()
