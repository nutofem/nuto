set(NuToMath_SWIG_MODULES
    NuToMath
    ModuleTestInterface
)

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

include_directories(${CMAKE_SOURCE_DIR}/src)

# build python modules
foreach(NuToMath_SWIG_MODULE ${NuToMath_SWIG_MODULES})
    nuto_swig_module(${NuToMath_SWIG_MODULE} ${NuToMath_SWIG_MODULE}.i Math)
endforeach()
