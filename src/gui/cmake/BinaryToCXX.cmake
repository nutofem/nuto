include(FindPythonInterp)
if(NOT PYTHONINTERP_FOUND)
  message(FATAL_ERROR "Python interpreter required")
endif()

# binary_to_cxx(output_file INPUT input_file [IDENTIFIER identifier])
function(binary_to_cxx output_file)
  set(options ${ARGN})

  list(GET options 0 keyword)
  if("${keyword}" STREQUAL "INPUT")
    list(GET options 1 input_name)
    list(REMOVE_AT options 0 1)
  else()
    message(FATAL_ERROR "expected INPUT parameter")
  endif()

  set(bin2cxx_args "${input_name}" "-o" "${output_file}")

  if(options)
    list(GET options 0 keyword)
    if("${keyword}" STREQUAL "IDENTIFIER")
      list(GET options 1 identifier)
      list(REMOVE_AT options 0 1)
      list(APPEND bin2cxx_args "-i" "${identifier}")
    endif()
  endif()

  find_file(bin2cxx_path "bin2cxx.py" PATHS ${CMAKE_MODULE_PATH}
	    NO_DEFAULT_PATH)

  add_custom_command(OUTPUT ${output_file}
		     COMMAND "${PYTHON_EXECUTABLE}" "${bin2cxx_path}" ${bin2cxx_args}
		     DEPENDS ${input_name})
  set_source_files_properties(${output_file} PROPERTIES GENERATED TRUE)
endfunction(binary_to_cxx)
