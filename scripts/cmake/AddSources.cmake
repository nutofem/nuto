function(add_sources)
    get_property(is_defined GLOBAL PROPERTY SourcesList DEFINED)
    if(NOT is_defined)
        define_property(GLOBAL PROPERTY SourcesList
            BRIEF_DOCS "List of source files"
            FULL_DOCS "List of source files to be compiled in one library")
    endif()

    # make absolute paths
    set(SRCS)
    foreach(s IN LISTS ARGN)
        if(NOT IS_ABSOLUTE "${s}")
            get_filename_component(s "${s}" ABSOLUTE)
        endif()
    list(APPEND SRCS "${s}")
    endforeach()

    # append to global list
    set_property(GLOBAL APPEND PROPERTY SourcesList "${SRCS}")
endfunction()
