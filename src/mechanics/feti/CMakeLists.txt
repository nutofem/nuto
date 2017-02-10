

add_library(NuToFeti FetiSolver.cpp StructureFeti.cpp)
install(TARGETS NuToFeti
        ARCHIVE DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH}
        LIBRARY DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH}
        )