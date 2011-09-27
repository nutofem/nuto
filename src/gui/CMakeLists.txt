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
find_package(wxWidgets COMPONENTS base core aui)
find_package_message(WX "Found wxWidgets: ${wxWidgets_FOUND}" ${wxWidgets_FOUND})
if (wxWidgets_FOUND)
  include(${wxWidgets_USE_FILE})
endif()

include(FindVTK)
if(VTK_FOUND)
  find_package_message(VTK "Found VTK: ${VTK_USE_FILE}" ${VTK_USE_FILE})
  include(${VTK_USE_FILE})
  try_compile(VTK_RECENT
	      ${CMAKE_CURRENT_BINARY_DIR}/
	      ${CMAKE_CURRENT_SOURCE_DIR}/cmake/vtk_recent.cpp
	      CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${VTK_INCLUDE_DIRS}" 
	      )
  if(NOT VTK_RECENT)
    set(VTK_FOUND FALSE)
  endif(NOT VTK_RECENT)
endif()

# Also need Boost.Python for GUI
if(Boost_FOUND)
  find_package(Boost ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION} COMPONENTS python)
endif()

if(wxWidgets_FOUND AND VTK_FOUND AND Boost_PYTHON_FOUND AND ENABLE_VISUALIZE)
  message(STATUS "All requirements for GUI met.")
  message(STATUS "\tFound wxWidgets: ${wxWidgets_FOUND}")
  message(STATUS "\t\t wxWidgets root dir    : ${wxWidgets_ROOT_DIR}")
  message(STATUS "\t\t wxWidgets include dirs: ${wxWidgets_INCLUDE_DIRS}")
  message(STATUS "\t\t wxWidgets library dirs: ${wxWidgets_LIBRARY_DIRS}")
  message(STATUS "\t\t wxWidgets libraries   : ${wxWidgets_LIBRARIES}")
  # Utilities
  include(ArtProvider)
  include(SplitDebugInfo)

  add_subdirectory(vtkwx)

  include_directories(include)
  include_directories(vtkwx)
  add_subdirectory(src)
else()
  message(STATUS "GUI not built: need: external libraries wxWidgets, VTK and Boost.Python")
  message(STATUS "... also, ENABLE_VISUALIZE must be set")
endif()
