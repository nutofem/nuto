add_subdirectory(base)
add_subdirectory(math)
if(ENABLE_MECHANICS)
    add_subdirectory(mechanics)
endif()

if(ENABLE_METAMODEL)
    add_subdirectory(metamodel)
endif()

if(ENABLE_OPTIMIZE)
    add_subdirectory(optimize)
endif()
