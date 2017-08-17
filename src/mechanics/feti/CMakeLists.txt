add_library(FetiLib FetiSolver.cpp StructureFeti.cpp NewmarkFeti.cpp)
set_target_properties(FetiLib PROPERTIES
    OUTPUT_NAME NuToFeti
    )
target_include_directories(FetiLib PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:include/nuto>
    )

install(TARGETS FetiLib EXPORT NuToTargets
    LIBRARY DESTINATION lib
    INCLUDES DESTINATION include/nuto
    )
