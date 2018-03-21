add_subdirectory(base)
add_subdirectory(math)
add_subdirectory(mechanics)
add_subdirectory(geometryConcrete)
add_subdirectory(visualize)

get_property(Sources GLOBAL PROPERTY SourcesList)
add_library(NuTo ${Sources})
target_link_libraries(NuTo PUBLIC Eigen3::Eigen)
target_link_libraries(NuTo PRIVATE ${LAPACK_LIBRARIES})
target_link_libraries(NuTo PRIVATE Boost::filesystem)
target_link_libraries(NuTo PRIVATE Rang)
if(SUITESPARSE_FOUND)
    target_link_libraries(NuTo PRIVATE SuiteSparse::UmfPack SuiteSparse::Cholmod)
endif()
if(MUMPS_FOUND)
    target_link_libraries(NuTo PRIVATE EigenMumpsSupport)
endif()
target_include_directories(NuTo PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/>
    $<INSTALL_INTERFACE:include/>
    )

install(TARGETS NuTo EXPORT NuToTargets
        LIBRARY DESTINATION lib
)
