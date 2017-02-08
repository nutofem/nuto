

add_library(NuToFeti FetiSolver.cpp NewmarkFeti.cpp StructureFETI.cpp)
install(TARGETS NuToFeti
        ARCHIVE DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH}
        LIBRARY DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH}
        )