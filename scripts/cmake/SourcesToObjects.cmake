# give it a list of sources, and it creates an object lib for each source,
# and returns the list of link targets
function(sources_to_objects Sources Objects)
    foreach(source ${Sources})
        create_object_lib(${CMAKE_CURRENT_SOURCE_DIR}/${source} object)
        list(APPEND objectList ${object})
    endforeach()
    set(${Objects} ${objectList} PARENT_SCOPE)
endfunction()
