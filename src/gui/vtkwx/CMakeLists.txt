# Check for GTK (needed when wxGTK is in use)
include(FindPkgConfig)
pkg_check_modules(GTK gtk+-2.0)

include_directories(${GTK_INCLUDE_DIRS})

try_compile(HAS_gtk_widget_get_window
            ${CMAKE_CURRENT_BINARY_DIR}/config.out 
	    ${CMAKE_CURRENT_SOURCE_DIR}/config/gtk_widget_get_window.c
	    CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${GTK_INCLUDE_DIRS}" "-DLINK_LIBRARIES=${GTK_LIBRARIES}")
find_package_message(vtkwx_HAS_gtk_widget_get_window
		     "gtk_widget_get_window() available: ${HAS_gtk_widget_get_window}"
		     ${HAS_gtk_widget_get_window})

if(HAS_gtk_widget_get_window)
  add_definitions(-DHAS_gtk_widget_get_window)
endif(HAS_gtk_widget_get_window)

# Library
file(GLOB vtkwx_sources *.h *.cpp)
add_library(vtkwx STATIC ${vtkwx_sources})
target_link_libraries(vtkwx vtkCommon vtkRendering)
# Link against GTK because we use gtk_widget_get_window() together with wxGTK
target_link_libraries(vtkwx ${GTK_LIBRARIES})
if(STATIC_LIBRARY_FLAGS)
  set_target_properties(vtkwx PROPERTIES COMPILE_FLAGS ${STATIC_LIBRARY_FLAGS})
endif()
