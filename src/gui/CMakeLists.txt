PROJECT(nutogui)

cmake_minimum_required(VERSION 2.8)

list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(CheckCXXSourceCompiles)
include(TestCXXAcceptsFlag)
# Basic checks
CHECK_CXX_SOURCE_COMPILES("int main() { char* p = nullptr; }" HAVE_NULLPTR)
if(HAVE_NULLPTR)
  add_definitions(-DHAVE_NULLPTR)
endif()

function(add_linker_flag flagstr FLAG)
  CHECK_CXX_ACCEPTS_FLAG(${flagstr} FLAG_LINKER_${FLAG})
  if(FLAG_LINKER_${FLAG})
    link_libraries(${flagstr})
  endif()
endfunction()
function(add_cxx_flag flagstr FLAG)
  CHECK_CXX_ACCEPTS_FLAG(${flagstr} FLAG_CXX_${FLAG})
  if(FLAG_CXX_${FLAG})
    add_definitions(${flagstr})
  endif()
endfunction()
add_linker_flag("-Wl,--as-needed" AS_NEEDED)
add_linker_flag("-Wl,--warn-unresolved-symbols" WARN_UNRESOLVED)
add_cxx_flag("-Wall" WARNINGS_ALL)
add_cxx_flag("-fvisibility-inlines-hidden" INLINES_VISIBILITY_HIDDEN)
add_cxx_flag("-fvisibility=hidden" DEFAULT_VISIBILITY_HIDDEN)

# These flags can make the final binary smaller
add_cxx_flag("-ffunction-sections" FUNCTION_SECTIONS)
add_cxx_flag("-fdata-sections" DATA_SECTIONS)
add_linker_flag("-Wl,--gc-sections" GC_SECTIONS)

if(BUILD_SHARED_LIBS)
  set(STATIC_LIBRARY_FLAGS -fPIC)
  add_definitions(-DNUTOGUI_SHARED_LIBS)
endif()

CHECK_CXX_ACCEPTS_FLAG("-Wl,--whole-archive -Wl,--no-whole-archive" FLAG_LINKER_WHOLE_ARCHIVE)
if(FLAG_LINKER_WHOLE_ARCHIVE)
  set(LINK_WHOLE_ARCHIVE "-Wl,--whole-archive")
  set(LINK_NO_WHOLE_ARCHIVE "-Wl,--no-whole-archive")
endif()

# Dependencies
include(FindwxWidgets)
find_package(wxWidgets REQUIRED base core aui)
include(${wxWidgets_USE_FILE})

include(FindVTK)
if(NOT VTK_FOUND)
  message(FATAL_ERROR "VTK is required.")
endif(NOT VTK_FOUND)
find_package_message(VTK "Found VTK: ${VTK_USE_FILE}" ${VTK_USE_FILE})
include(${VTK_USE_FILE})
try_compile(VTK_RECENT
            ${CMAKE_CURRENT_BINARY_DIR}/
            ${CMAKE_CURRENT_SOURCE_DIR}/cmake/vtk_recent.cpp
            CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${VTK_INCLUDE_DIRS}" 
            )
if(NOT VTK_RECENT)
  message(FATAL_ERROR "VTK is too old. Version 5.2 required.")
endif(NOT VTK_RECENT)

# Boost: as of 2010-11-03, NuTo requires Boost 1.44.0; hence, so do we...
include(FindBoost)
SET(Boost_ADDITIONAL_VERSIONS)
find_package(Boost 1.44.0 REQUIRED python serialization)
include_directories(${Boost_INCLUDE_DIRS})

# @@@ Write a proper check ...
if(NOT NUTO_PREFIX)
  message(FATAL_ERROR "NUTO_PREFIX variable not specified.")
endif()
get_filename_component(NUTO_PREFIX "${NUTO_PREFIX}" ABSOLUTE)
find_package_message(NUTO "Using NUTO: ${NUTO_PREFIX}" ${NUTO_PREFIX})
include_directories(${NUTO_PREFIX}/include)
link_libraries(-L${NUTO_PREFIX}/lib)

# Utilities
include(ArtProvider)
include(SplitDebugInfo)

add_subdirectory(vtkwx)

include_directories(include)
include_directories(vtkwx)
add_subdirectory(src)
