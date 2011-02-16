#include(CMakeFindBinUtils)
find_program(CHMOD chmod)

# split_debug_info - uses objcopy to separate debug info of a target into a separate file
function (split_debug_info target)
  if(CMAKE_OBJCOPY)
    get_target_property(${target}_LOCATION ${target} LOCATION)
    set(${target}_DEBUG_FILE "${${target}_LOCATION}.dbg")
    add_custom_command(TARGET ${target}
		      POST_BUILD
		      COMMAND ${CMAKE_OBJCOPY} ARGS --only-keep-debug "${${target}_LOCATION}" "${${target}_DEBUG_FILE}"
		      COMMAND ${CMAKE_OBJCOPY} ARGS --strip-unneeded "${${target}_LOCATION}"
		      COMMAND ${CMAKE_OBJCOPY} ARGS "--add-gnu-debuglink=${${target}_DEBUG_FILE}" "${${target}_LOCATION}" # Need --long-section-names argument on some platforms
		      COMMAND ${CHMOD} ARGS a-x "${${target}_DEBUG_FILE}"
		      COMMENT "Separating debug info from ${${target}_LOCATION}")
  endif()
endfunction(split_debug_info)
