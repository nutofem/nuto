add_library(NuToFeti FetiSolver.cpp StructureFeti.cpp NewmarkFeti.cpp)

install(TARGETS NuToFeti
    LIBRARY DESTINATION lib
    )
