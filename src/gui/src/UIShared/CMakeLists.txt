add_library(UIShared dummy.cpp)

# Link other implementation libraries into UIShared.
# Done in a roundabout way so the .a files are enclosed in 'whole-archive'
# flags, but none of the dependencies.
set(UIShared_deps Platform UICommon GuiFrame ResultViewer ScriptEditor)
# LIB_DEPENDS is a list, each two entries make a pair.
# The first entry sets the "link type" for an argument (2nd entry).
# The link type can be general, optimized or debug.
# The code below convertes the CMAKE_BUILD_TYPE to that link type.
get_property(DEBUG_CONFIGURATIONS GLOBAL PROPERTY DEBUG_CONFIGURATIONS)
if(NOT "${CMAKE_BUILD_TYPE}" STREQUAL "")
  string(TOUPPER "${CMAKE_BUILD_TYPE}" build_type_upper)
  if("${DEBUG_CONFIGURATIONS}")
    foreach(debug_config ${DEBUG_CONFIGURATIONS})
      string(TOUPPER "${debug_config}" debug_config_upper)
      if("${build_type_upper}" STREQUAL "${debug_config_upper}")
	set(effective_build_type debug)
	break()
      endif()
    endforeach()
  else()
    if("${build_type_upper}" STREQUAL "DEBUG")
      set(effective_build_type "debug")
    endif()
  endif()
endif()
if("${effective_build_type}" STREQUAL "")
  set(effective_build_type "optimized")
endif()
# Actually link libraries.
foreach(dep ${UIShared_deps})
  # Actual library
  get_target_property(${dep}_LOCATION ${dep} LOCATION)
  target_link_libraries(UIShared "${LINK_WHOLE_ARCHIVE} ${${dep}_LOCATION} ${LINK_NO_WHOLE_ARCHIVE}")
  # FIXME: Crude way to make UIShared depend on dep - using add_dependencies() didn't seem to have the desired effect...
  target_link_libraries(UIShared ${dep})
  # Dependencies
  set(dep_flags ${${dep}_LIB_DEPENDS})
  while(dep_flags)
    list(GET dep_flags 0 dep_flag_config)
    list(GET dep_flags 1 dep_flag_actual)
    list(REMOVE_AT dep_flags 0 1)
    if(("${dep_flag_config}" STREQUAL "${effective_build_type}")
	OR ("${dep_flag_config}" STREQUAL "general"))
      target_link_libraries(UIShared ${dep_flag_actual})
    endif()
  endwhile()
endforeach()

target_link_libraries(UIShared ${wxWidgets_LIBRARIES})
set_target_properties(UIShared
		      PROPERTIES
		      LINK_INTERFACE_LIBRARIES ""
		      LIBRARY_OUTPUT_NAME "nutogui_UIShared")
if(BUILD_SHARED_LIBS)
  install(TARGETS UIShared
	  RUNTIME DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH}
	  LIBRARY DESTINATION ${NUTO_LIBRARIES_INSTALL_PATH})
endif()
