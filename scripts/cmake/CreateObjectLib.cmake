function(create_object_lib Source Object)
    # give a source file, and create an object lib for it
    # return the link target $<TARGET_OBJECT:target> to Object
    get_target_from_source(${Source} libName)
    add_library(${libName} OBJECT ${Source})
    target_include_directories(${libName} PUBLIC ${CMAKE_SOURCE_DIR}/src)

    set(objectLib "$<TARGET_OBJECTS:${libName}>")
    set(${Object} ${objectLib} PARENT_SCOPE)
endfunction()
