# temporarily save compiler flags, set them to "python flags", and reset them
# at the end of the file
set(ORIGINAL_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "${PYTHON_CXX_FLAGS}")

find_package(Numpy REQUIRED)
include_directories(${NUMPY_INCLUDE_PATH})

# Python required for interface files c++ and python
message(STATUS "Checking for Python...")
find_package(PythonInterp 3 REQUIRED)
find_package(PythonLibs 3 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
# get python version
execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import sys;print(sys.version[:3]);"
        OUTPUT_VARIABLE PYTHON_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS "PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
message(STATUS "PYTHON_INCLUDE_DIRS = ${PYTHON_INCLUDE_DIRS}")

# SWIG required for interface files c++ and python
message(STATUS "Checking for SWIG...")
find_program(SWIG_EXECUTABLE NAMES swig3.0 swig swig2.0)
message(STATUS "SWIG_EXECUTABLE = ${SWIG_EXECUTABLE1}")
find_package(SWIG REQUIRED)
message(STATUS "SWIG found: ${SWIG_EXECUTABLE}")
message(STATUS "SWIG_VERSION = ${SWIG_VERSION}")
if(SWIG_VERSION VERSION_LESS "3.0")
  message(FATAL_ERROR "Insufficient swig version (at least 3.0)")
endif()
include(${SWIG_USE_FILE})
set(NuTo_SWIG_FLAGS "-Wall;-fvirtual;${NuTo_SWIG_FLAGS}")

configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py)
include_directories(${CMAKE_SOURCE_DIR}/src)

add_subdirectory(base)
add_subdirectory(math)
if(ENABLE_MECHANICS)
    add_subdirectory(mechanics)
else()
    deactivate_module(mechanics.ModuleMechanics)
endif()

if(ENABLE_GEOMETRYCONCRETE)
    add_subdirectory(geometryConcrete)
else()
    deactivate_module(geometryConcrete.ModuleGeometryConcrete)
endif()

if(ENABLE_METAMODEL)
    add_subdirectory(metamodel)
else()
    deactivate_module(metamodel.ModuleMetamodel)
endif()

if(ENABLE_OPTIMIZE)
    add_subdirectory(optimize)
else()
    deactivate_module(optimize.ModuleOptimize)
endif()

set(CMAKE_CXX_FLAGS "${ORIGINAL_CXX_FLAGS}")
