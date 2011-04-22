add_subdirectory(ScriptRunner)

file(GLOB app_art art/tango-icons/*.png)
bundle_artprovider_art(app_art_sources "image/png" ${app_art})

set(THROBBER_SOURCE ${CMAKE_CURRENT_BINARY_DIR}/throbber.png.inc)
binary_to_cxx(${THROBBER_SOURCE}
	      INPUT ${CMAKE_CURRENT_SOURCE_DIR}/art/tango-process-working.png
	      IDENTIFIER throbber_png)
set_source_files_properties(${THROBBER_SOURCE} PROPERTIES GENERATED TRUE)

file(GLOB app_sources *.h *.cpp)
include_directories(${CMAKE_CURRENT_BINARY_DIR})
ADD_EXECUTABLE(nutogui ${app_sources} ${app_art_sources} ${THROBBER_SOURCE})
set_target_properties(nutogui
		      PROPERTIES
		      RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
target_link_libraries(nutogui UIShared ScriptRunner)
target_link_libraries(nutogui vtkIO) # for vtkDataSetReader
split_debug_info(nutogui)
install(TARGETS nutogui DESTINATION bin)
