# Check if wxStyledTextCtrl is available.
# It's core on WX 2.9, but not 2.8, where it's a 'contrib' module and not generally available.
# In that case, link our own copy (from WX 2.8.11).
execute_process(COMMAND ${wxWidgets_CONFIG_EXECUTABLE} --libs stc
		OUTPUT_VARIABLE wxSTC_CONFIG_LIBS
		OUTPUT_STRIP_TRAILING_WHITESPACE)

# try_compile's COMPILE_DEFINITIONS options wants -Dfoo, but wxWidgets_DEFINITIONS is bare
foreach(wxdef ${wxWidgets_DEFINITIONS})
  list(APPEND wxWidgets_DEFINITIONS_PROCESSED -D${wxdef})
endforeach()

try_compile(HAS_WX_STC
            ${CMAKE_CURRENT_BINARY_DIR}/
	    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/wxstc.cpp
	    CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${wxWidgets_INCLUDE_DIRS};${CMAKE_CURRENT_SOURCE_DIR}/wxstc/include" 
	      "-DLINK_LIBRARIES=${wxWidgets_LIBRARIES};${wxSTC_CONFIG_LIBS}"
	      "-DLINK_DIRECTORIES=${wxWidgets_LIBRARY_DIRS}"
	    COMPILE_DEFINITIONS ${wxWidgets_DEFINITIONS_PROCESSED}
	    )
find_package_message(HAS_WX_STC
		     "WX has wxStyledTextControl: ${HAS_WX_STC}"
		     ${HAS_WX_STC})
if(HAS_WX_STC)
  # Use existing wxSTC
  set(wxSTC_LIBRARIES ${wxSTC_CONFIG_LIBS})
else()
  # Build our own copy of wxSTC
  add_subdirectory(wxstc)
  set(wxSTC_LIBRARIES wxstc)
  include_directories("${CMAKE_CURRENT_SOURCE_DIR}/wxstc/include")
endif()

add_definitions(-DSCRIPTEDITOR_EXPORTED=NUTOGUI_EXPORT)
file(GLOB ScriptEditor_sources *.h *.cpp)
add_library(ScriptEditor STATIC ${ScriptEditor_sources})
target_link_libraries(ScriptEditor ${wxSTC_LIBRARIES})
target_link_libraries(ScriptEditor UICommon)
if(STATIC_LIBRARY_FLAGS)
  set_target_properties(ScriptEditor PROPERTIES COMPILE_FLAGS ${STATIC_LIBRARY_FLAGS})
endif()
