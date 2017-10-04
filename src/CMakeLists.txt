add_subdirectory(base)
add_subdirectory(math)

add_subdirectory(mechanics)

if(ENABLE_GEOMETRYCONCRETE)
    add_subdirectory(geometryConcrete)
endif()

add_subdirectory(metamodel)

add_subdirectory(optimize)

if(ENABLE_VISUALIZE)
    add_subdirectory(visualize)
endif()
