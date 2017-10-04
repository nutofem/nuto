# temporarily save compiler flags, set them to "python flags", and reset them
# at the end of the file
set(ORIGINAL_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "${PYTHON_CXX_FLAGS}")

# Python required for interface files c++ and python
message(STATUS "Checking for Python...")
find_package(PythonInterp 3 REQUIRED)
find_package(PythonLibs 3 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
# get python version
execute_process(
    COMMAND ${PYTHON_EXECUTABLE} -c "import sys;print(sys.version[:3]);"
    OUTPUT_VARIABLE PYTHON_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS "PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
message(STATUS "PYTHON_INCLUDE_DIRS = ${PYTHON_INCLUDE_DIRS}")
message(STATUS "PYTHON_VERSION = ${PYTHON_VERSION}")

# set numpy include path for generated wrapper files
execute_process(COMMAND
    ${CMAKE_SOURCE_DIR}/scripts/get_numpy_path.py
    OUTPUT_VARIABLE NUMPY_PATH
    OUTPUT_STRIP_TRAILING_WHITESPACE)
set(PYTHON_CXX_FLAGS "${PYTHON_CXX_FLAGS} -I${NUMPY_PATH}")

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
set(NuTo_SWIG_FLAGS "-fvirtual;${NuTo_SWIG_FLAGS}")

configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py)
include_directories(${CMAKE_SOURCE_DIR}/src)

add_subdirectory(base)
add_subdirectory(math)
add_subdirectory(mechanics)

if(ENABLE_GEOMETRYCONCRETE)
    add_subdirectory(geometryConcrete)
else()
    deactivate_module(geometryConcrete.ModuleGeometryConcrete)
endif()

add_subdirectory(metamodel)

add_subdirectory(optimize)

# support for strongly typed enums was introduced in 3.0.3
if(SWIG_VERSION VERSION_LESS "3.0.3")
    file(READ ${CMAKE_CURRENT_BINARY_DIR}/__init__.py FILE_CONTENT)
    file(READ enum_hack.py ADDITIONAL_ENUMS)
    string(CONCAT MODIFIED_FILE_CONTENT ${FILE_CONTENT} ${ADDITIONAL_ENUMS})
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/__init__.py
        "${MODIFIED_FILE_CONTENT}")
endif()

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    DESTINATION lib/python${PYTHON_VERSION}/dist-packages
    FILES_MATCHING PATTERN "*.py"
    PATTERN "*.so"
    PATTERN "*.dir" EXCLUDE
    PATTERN "CMakeFiles" EXCLUDE
    )

set(CMAKE_CXX_FLAGS "${ORIGINAL_CXX_FLAGS}")
