set(NuToOptimizeSources
    ConjugateGradientGrid.cpp
    ConjugateGradientNonLinear.cpp
    CallbackHandler.cpp
    CallbackHandlerPython.cpp
    CallbackHandlerGrid.cpp
    Jacobi.cpp
    MisesWielandt.cpp
    NonlinearSolverBase.cpp
    NewtonRaphson.cpp
    )

message(STATUS "Checking for Python...")
find_package(PythonInterp 3 REQUIRED)
find_package(PythonLibs 3 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

sources_to_objects("${NuToOptimizeSources}" NuToOptimizeObjects)

add_library(NuToOptimize ${NuToOptimizeSources})
target_link_libraries(NuToOptimize NuToMath NuToBase ${PYTHON_LIBRARY})

if(ENABLE_VISUALIZE)
    target_link_libraries(NuToOptimize NuToVisualize)
endif()

if(ENABLE_SHOW_TIME)
    #add library for clock_gettime routines
    target_link_libraries(NuToOptimize rt)
endif()

set_target_properties(NuToOptimize PROPERTIES COMPILE_FLAGS "${PYTHON_C_FLAGS}")

# install library
install(TARGETS NuToOptimize
    ARCHIVE DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH}
    LIBRARY DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH}
    )
