find_package(Numpy REQUIRED)
include_directories(${NUMPY_INCLUDE_PATH})

# Python required for interface files c++ and python
message(STATUS "Checking for Python...")
find_package(PythonInterp 2 REQUIRED)
find_package(PythonLibs 2 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
# get python version
execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import sys;print(sys.version[:3]);"
        OUTPUT_VARIABLE PYTHON_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS "PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
message(STATUS "PYTHON_INCLUDE_DIRS = ${PYTHON_INCLUDE_DIRS}")

set(NUTO_PYTHON_MODULES_INSTALL_PATH "lib${LIB_SUFFIX}/python${PYTHON_VERSION}/site-packages")

# SWIG required for interface files c++ and python
message(STATUS "Checking for SWIG...")
find_program(SWIG_EXECUTABLE NAMES swig3.0 swig swig2.0)
message(STATUS "MYSWIG_EXECUTABLE = ${SWIG_EXECUTABLE1}")
find_package(SWIG REQUIRED)
message(STATUS "SWIG found: ${SWIG_EXECUTABLE}")
message(STATUS "SWIG_VERSION = ${SWIG_VERSION}")
if(SWIG_VERSION VERSION_LESS "3.0")
  message(FATAL_ERROR "Insufficient swig version (at least 3.0)")
endif()
include(${SWIG_USE_FILE})
set(NuTo_SWIG_FLAGS "-Wall;-fvirtual;${NuTo_SWIG_FLAGS}")

configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py)
nuto_install_python_file(__init__.py ${NUTO_PYTHON_MODULES_INSTALL_PATH}/nuto)
include_directories(${CMAKE_SOURCE_DIR}/src)

add_subdirectory(base)
add_subdirectory(math)
if(ENABLE_MECHANICS)
    add_subdirectory(mechanics)
else()
    set(SEARCH_REGEX "from nuto.mechanics.ModulMechanics")
    file(READ ${CMAKE_CURRENT_BINARY_DIR}/__init__.py FILE_CONTENT)
    string(REGEX REPLACE "${SEARCH_REGEX}" "#from nuto.mechanics.ModulMechanics"
        MODIFIED_FILE_CONTENT "${FILE_CONTENT}")
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/__init__.py "${MODIFIED_FILE_CONTENT}")
endif()

if(ENABLE_GEOMETRYCONCRETE)
    add_subdirectory(geometryConcrete)
else()
    set(SEARCH_REGEX "from nuto.geometryConcrete.ModulGeometryConcrete")
    file(READ ${CMAKE_CURRENT_BINARY_DIR}/__init__.py FILE_CONTENT)
    string(REGEX REPLACE "${SEARCH_REGEX}" "#from nuto.geometryConcrete.ModulGeometryConcrete"
        MODIFIED_FILE_CONTENT "${FILE_CONTENT}")
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/__init__.py "${MODIFIED_FILE_CONTENT}")
endif()

if(ENABLE_METAMODEL)
    add_subdirectory(metamodel)
else()
    set(SEARCH_REGEX "from nuto.metamodel.ModulMetaModel")
    file(READ ${CMAKE_CURRENT_BINARY_DIR}/__init__.py FILE_CONTENT)
    string(REGEX REPLACE "${SEARCH_REGEX}" "#from nuto.metamodel.ModulMetaModel"
        MODIFIED_FILE_CONTENT "${FILE_CONTENT}")
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/__init__.py "${MODIFIED_FILE_CONTENT}")
endif()

if(ENABLE_OPTIMIZE)
    add_subdirectory(optimize)
else()
    set(SEARCH_REGEX "from nuto.optimize.ModulOptimizer")
    file(READ ${CMAKE_CURRENT_BINARY_DIR}/__init__.py FILE_CONTENT)
    string(REGEX REPLACE "${SEARCH_REGEX}" "#from nuto.optimize.ModulOptimizer"
        MODIFIED_FILE_CONTENT "${FILE_CONTENT}")
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/__init__.py "${MODIFIED_FILE_CONTENT}")
endif()
