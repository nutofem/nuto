add_subdirectory(base)
add_subdirectory(math)

add_subdirectory(mechanics)

if(ENABLE_GEOMETRYCONCRETE)
    add_subdirectory(geometryConcrete)
endif()

if(ENABLE_METAMODEL)
    add_subdirectory(metamodel)
endif()

if(ENABLE_OPTIMIZE)
    add_subdirectory(optimize)
endif()

if(ENABLE_VISUALIZE)
    add_subdirectory(visualize)
endif()
