# nuto_swig_module(<module> <interface file> <install_path> <dependencies>)
#   Function to set up a SWIG module for a NuTo library.
#   <module>        SWIG module name
#   <interface file>    SWIG .i file
#   <dependencies>      Dependencies needed by that module.
#                       Typically the NuTo library it wraps.
function(nuto_swig_module module_name interface_file)
    set(libraries ${ARGN})

    include_directories(SYSTEM ${CMAKE_SOURCE_DIR}/external)
    set_source_files_properties(${interface_file} PROPERTIES CPLUSPLUS ON)
    set_source_files_properties(${interface_file}
        PROPERTIES SWIG_FLAGS "${NuTo_SWIG_FLAGS}")
    swig_add_module(${module_name} python ${interface_file})
    # link library
    swig_link_libraries(${module_name} ${libraries} ${PYTHON_LIBRARIES})
    # check for unresolved symbols
    set_target_properties(${SWIG_MODULE_${module_name}_REAL_NAME}
        PROPERTIES LINK_FLAGS -Wl,-z,defs)
    # Additional build flags for module
    set_source_files_properties("${swig_generated_file_fullname}"
        PROPERTIES COMPILE_FLAGS "${PYTHON_CXX_FLAGS}")
endfunction()

