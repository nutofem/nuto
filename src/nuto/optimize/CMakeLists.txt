set(NuToOptimize_SOURCES
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

# build shared library
add_library(NuToOptimize ${NuToOptimize_SOURCES})
target_link_libraries(NuToOptimize NuToMath NuToBase ${PYTHON_LIBRARY})

if(ENABLE_VISUALIZE)
    target_link_libraries(NuToOptimize NuToVisualize)
endif()

if(ENABLE_SHOW_TIME)
    #add library for clock_gettime routines
    target_link_libraries(NuToOptimize rt)
endif()

set(NuToOptimize_COMPILE_FLAGS "${PYTHON_C_FLAGS}")
if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
    set(NuToOptimize_COMPILE_FLAGS "${NuToOptimize_COMPILE_FLAGS} -fPIC")
endif()
set_target_properties(NuToOptimize PROPERTIES COMPILE_FLAGS "${NuToOptimize_COMPILE_FLAGS}")

# install library
install(TARGETS NuToOptimize
    ARCHIVE DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH}
    LIBRARY DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH}
    )

# build python modules
configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py)
nuto_install_python_file(__init__.py ${NUTO_PYTHON_MODULES_INSTALL_PATH}/nuto/optimize)
nuto_swig_module(ModulOptimizer ModulOptimizer.i nuto/optimize NuToOptimize)
