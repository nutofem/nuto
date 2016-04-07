SET(NuToOptimize_SOURCES
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

SET(NuToOptimize_SWIG_MODULES
  ModulOptimizer
)

# build shared library
ADD_LIBRARY(NuToOptimize ${NuToOptimize_SOURCES})
TARGET_LINK_LIBRARIES(NuToOptimize NuToMath NuToBase)
IF (ENABLE_MECHANICS)
  TARGET_LINK_LIBRARIES(NuToOptimize NuToMechanics)
ENDIF()
IF (ENABLE_VISUALIZE)
  TARGET_LINK_LIBRARIES(NuToOptimize NuToVisualize)
ENDIF (ENABLE_VISUALIZE)
IF (ENABLE_SHOW_TIME)
#add library for clock_gettime routines
  TARGET_LINK_LIBRARIES(NuToOptimize rt)
ENDIF (ENABLE_SHOW_TIME)

set(NuToOptimize_COMPILE_FLAGS "${PYTHON_C_FLAGS}")
IF(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
  set(NuToOptimize_COMPILE_FLAGS "${NuToOptimize_COMPILE_FLAGS} -fPIC")
ENDIF(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
SET_TARGET_PROPERTIES(NuToOptimize PROPERTIES COMPILE_FLAGS "${NuToOptimize_COMPILE_FLAGS}")

# install library
INSTALL(TARGETS NuToOptimize
	ARCHIVE DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH}
	LIBRARY DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH})

# build python modules
CONFIGURE_FILE(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py)
NUTO_INSTALL_PYTHON_FILE(__init__.py ${NUTO_PYTHON_MODULES_INSTALL_PATH}/nuto/optimize)
FOREACH(NuToOptimize_SWIG_MODULE ${NuToOptimize_SWIG_MODULES})
  NUTO_SWIG_MODULE(${NuToOptimize_SWIG_MODULE} ${NuToOptimize_SWIG_MODULE}.i nuto/optimize NuToOptimize)
ENDFOREACH(NuToOptimize_SWIG_MODULE)
