set(OptimizeSources
    ConjugateGradientNonLinear.cpp
    CallbackHandler.cpp
    CallbackHandlerPython.cpp
    Jacobi.cpp
    MisesWielandt.cpp
    NonlinearSolverBase.cpp
    NewtonRaphson.cpp
    )

message(STATUS "Checking for Python...")
find_package(PythonInterp 3 REQUIRED)
find_package(PythonLibs 3 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

sources_to_objects("${OptimizeSources}" OptimizeObjects)

add_library(Optimize ${OptimizeSources})
set_target_properties(Optimize PROPERTIES
    OUTPUT_NAME NuToOptimize
    )

target_link_libraries(Optimize Math Base ${PYTHON_LIBRARY})

if(ENABLE_VISUALIZE)
    target_link_libraries(Optimize Visualize)
endif()

set_target_properties(Optimize PROPERTIES COMPILE_FLAGS "${PYTHON_C_FLAGS}")

install(TARGETS Optimize EXPORT NuToTargets
    LIBRARY DESTINATION lib
    )
