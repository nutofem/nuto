if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    add_definitions(-Wsuggest-override)
endif()


add_subdirectory(base)
add_subdirectory(math)

if(ENABLE_MECHANICS)
    add_subdirectory(mechanics)
endif()

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
