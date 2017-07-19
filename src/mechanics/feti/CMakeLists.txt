add_library(FetiLib FetiSolver.cpp StructureFeti.cpp NewmarkFeti.cpp)
set_target_properties(FetiLib PROPERTIES
    OUTPUT_NAME NuToFeti
    )

install(TARGETS FetiLib EXPORT NuToTargets
    LIBRARY DESTINATION lib
    )
