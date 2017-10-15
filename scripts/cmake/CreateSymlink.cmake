# symlink a given file/directory to the build directory
function(create_symlink pathName)
    execute_process(COMMAND "${CMAKE_COMMAND}" "-E" "create_symlink"
        "${CMAKE_CURRENT_SOURCE_DIR}/${pathName}"
        "${CMAKE_CURRENT_BINARY_DIR}/${pathName}")
endfunction()
